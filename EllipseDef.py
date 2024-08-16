#!/usr/bin/python
# -*- coding: utf-8 -*-
#Solve for the best fit for an ellipse then plot it!
# Generated: 26 July 2024 E Fernandez

##This is a continuation off the fitEllipse3_new.py python file wherein the actual ellipse calculation definitions are included here. 

#Follows an approach suggested by Fitzgibbon, Pilu and Fischer in Fitzgibbon, A.W., Pilu, M., 
# and Fischer R.B., Direct least squares fitting of ellipsees, Proc. of the 13th Internation 
# Conference on Pattern Recognition, pp 253–257, Vienna, 1996.  
# Discussed on http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html and uses relationships 
#  found at http://mathworld.wolfram.com/Ellipse.html.
#
# ***Update from Authors: Andrew Fitzgibbon, Maurizio Pilu, Bob Fisher
# http://research.microsoft.com/en-us/um/people/awf/ellipse/
# Reference: "Direct Least Squares Fitting of Ellipses", IEEE T-PAMI, 1999
# Citation:  Andrew W. Fitzgibbon, Maurizio Pilu, and Robert B. Fisher
# Direct least-squares fitting of ellipses,
# IEEE Transactions on Pattern Analysis and Machine Intelligence, 21(5), 476--480, May 1999
#  @Article{Fitzgibbon99,
#   author = "Fitzgibbon, A.~W.and Pilu, M. and Fisher, R.~B.",
#   title = "Direct least-squares fitting of ellipses",
#   journal = pami,
#   year = 1999,
#   volume = 21,
#   number = 5,
#   month = may,
#   pages = "476--480"
#  }
# 
# This is a more bulletproof version than that in the paper, incorporating
# scaling to reduce roundoff error, correction of behaviour when the input 
# data are on a perfect hyperbola, and returns the geometric parameters
# of the ellipse, rather than the coefficients of the quadratic form.


import datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.feature
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs
import xarray as xr
import math
import netCDF4
#from get_ellipse_metrics import get_emetrics_max_min
from fitEllipse3_new import fitEllipseContour
from fitEllipse3_new import point_inside_polygon
#from EllipseCalc import ellipse_calc
from geopy.distance import great_circle
from matplotlib.patches import Polygon
from matplotlib import cm #colormaps!
import copy
import os
import scipy.stats
import pickle

##first two sets of definitions save for ERA5 data
##second two are for ERAI data

#in all loops there are commented out savefig statements for the ellipse images. 
#Please uncomment these and note your respective save location to use these images. 

########################################################################################################################
########################################################################################################################
##definition statement for calculating ellipse diagnostics with means/medians taken for multiple vortices

##between both definition statements the inputs may change in specification, and the empty lists for the variables must be noted ahead of time. 

def ERA5_one_ellipse_calc(lev_list,contour_list,windlat,temp_lev,y,m,t_d,
                     rat10,wind10,cenlat10,cenlon10,size10,ep10,numb10,
                     rat50,wind50,cenlat50,cenlon50,size50,ep50,numb50):
    #Start with defining desired period of record
    year = y
    month = m
    day = 1
    date1 = dt.datetime(year,month,day,0)  #first date to plot
    total_days = t_d
    first_fhr = 0
    hours = total_days * 24
    hr_inc = 24
    times = [date1 + dt.timedelta(hours=x) for x in range(0,hours,hr_inc)]
    date_list = netCDF4.date2num(times,units="hours since 1800-01-01 00:00:00",calendar="gregorian") #change dates to netcdf times

    ##start loop for calculation and plotting
    for plot_lev,the_contour in zip(lev_list,contour_list): #loop through each level
        print("Opening ERA5 data "+str(plot_lev)+" and contour level "+str(the_contour)  )

        #empty lists for appending diagnostics
        rat = []
        cenlt = []
        cenln = []
        wind = []
        sz = []
        ep = []
        num = []
        
        #text_file2 = "./"+date1.strftime("%Y%m%d_%H%M")+"_temp"+str(plot_lev)+".txt"
        #fout2 = open(text_file2,'w')
        if plot_lev == 10:
            ###THESE LINES MUST BE CHANGED FOR WHERE YOUR DATA IS STORED
            gfile = xr.open_dataset("../era5/gph/era5_gph_"+str(year)+".nc")
            g_files = gfile["z"]
            tfile = xr.open_dataset("../era5/t/era5_t_"+str(year)+".nc")
            t_files = tfile["t"]
            ufile = xr.open_dataset("../era5/u/era5_uv_"+str(year)+".nc")
            u_files = ufile["u"]
            
            print("Partitioning and averaging data.")
            print("First, g data")
            g_data = g_files.loc[dict(latitude=slice(90,0),level=lev_list)]
            print("u data")
            u_data = u_files.loc[dict(latitude=slice(windlat+1,windlat-1),level=lev_list)].mean(dim='longitude')
            print("t data")
            t_data = t_files.loc[dict(latitude=slice(90,75),level=lev_list)]
            tlats = t_data['latitude'].values
            weights=np.cos(np.deg2rad(t_data.latitude))
            weights.name = "weights"
            
            print("lat and lon")
            lats_era = g_data["latitude"].values
            lons_era = g_data["longitude"].values
        
        print(f"level: {plot_lev}")
        for date in date_list: #loop through each day 
            print(f"date: {date_list}")
            ind = np.where(date_list == date)[0]
            for d in ind:
                t = d
                height = g_data.loc[dict(time=times[t],level=plot_lev)] ##convert to hPa, gph
                height = height/9.81
                u = u_data.loc[dict(time=times[t],level=plot_lev)].mean(dim='latitude')
                tdata = t_data.loc[dict(time=times[t],level=plot_lev)].mean(dim='longitude')
            temp = tdata.weighted(weights).mean(dim='latitude')
            
            formatted_date = times[t].strftime("%Y%m%d")
            emark = []
            eline = []
            cs_temp = []
            print("##########################")
            print(date)
            print("Date "+times[t].strftime("%H UTC %d %b %Y"))
            valid_label = (times[t].strftime("Valid: %H UTC %d %b %Y"))
            
            plt.Figure(figsize=(15,15),dpi=120)#figure(figsize=(12,12),dpi=1200) <---Set fig size here!
            #ax = plt.subplot(1,1,1,projection=ccrs.Orthographic(0,90))
            ax = plt.axes(projection=ccrs.Orthographic(0,90))
            ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='dimgray',facecolor='none')
            #ax.set_extent([-180,180,0,90],ccrs.PlateCarree())
            ax.outline_patch.set_edgecolor('none')
            gl = ax.gridlines(color="grey",linestyle=":",linewidth=0.5)#=lat_lon_width)
            
            clevs = range(18000,33500,250)
            
            [x,y] = np.meshgrid(lons_era,lats_era)
            
            plt.title(str(plot_lev)+"hPa Elliptical Diagnostics\n"+valid_label)
            
            mem_color = "blue"
            mem_lw = 3.0
            mem_ms = 8
            cont_color = "#aec2e2"
            cont_lw = 1.0
            
            cs_temp.append(ax.contour(x,y,height,levels=clevs,linewidths=cont_lw,colors=cont_color))
            
            line_w = np.array(clevs,dtype=np.float)
            for i,c in enumerate(clevs):
                if c == the_contour:
                    line_w[i] = 1.5 # Make the contour we're using for the ellipse calc fatter
                else:
                    line_w[i] = 0.5
            print(clevs)
            ax.contour(x,y,height,clevs,transform=ccrs.PlateCarree(),extend='both',colors='black',linewidths=line_w)
            #check number/presence of contours
            try:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                print("in try")
                print("Number of ellipses: ",len(isoline_list))
            except:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []
            
            if len(isoline_list) > 0:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                #print(isoline_list)
                print("in ellipse ...")
                #quit()
                print("Number of ellipses: ",len(isoline_list))
            else:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []

            ##number of contours
            isocount = 0
            small = 0
            #empty lists for upcoming metric calculations
            ratio1 = []
            uvalues1 = []
            cenlat1 = []
            cenlon1 = []
            size1 = []
            ephi1 = []
            
            for isoline in isoline_list:
                ##next few lines will check if the contour crossing the prime meridian is one continuous contour or two separate ones
                #print('isoline',isoline)
                #[iso_lon,iso_lat] = mm(isoline[:,0],isoline[:,1],inverse=True)
                [iso_lon,iso_lat] = [isoline[:,0],isoline[:,1]]
                if len(iso_lon)<15:
                    print("-----Not analyzing ellipse with",len(iso_lon),"points, continuing...")  # Check for size!
                    small = 1
                    continue
                # Checking to see if contours are closed (0) or need to be joined (1)  Checking before convert lat/lon to radians
                lon_diff = abs(iso_lon[0] - iso_lon[len(iso_lon)-1])
                lat_diff = abs(iso_lat[0] - iso_lat[len(iso_lat)-1])
                join = 0
                if lon_diff > 1 or lat_diff > 1:
                    join = 1
                    print("Diffs lat/lon: ",lat_diff,lon_diff)
                iso_lon = np.deg2rad(iso_lon)
                iso_lat = np.deg2rad(iso_lat)
                
                ex = np.array((np.cos(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ey = np.array((np.sin(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ## does this contour include the pole?
                overpole = point_inside_polygon(0,0,ex,ey)  #returns true if poly includes the pole, false if not 
                print("It is",overpole, "that the contour includes the pole")
                if not overpole and len(isoline_list) > 1 and join > 0 and small < 1:
                    if isocount > 0:
                        ex = np.append(ex,ex2)
                        ey = np.append(ey,ey2)
                        print("Point didn't include pole - add to these to previous set")
                    else:
                        ex2 = ex
                        ey2 = ey
                        print("Point didn't include pole - keeping these to add to next set")
                        isocount = 1
                        continue
                ##begin actual calculation        
                print("Running ellipse diagnostic now")
                exx,eyy,eaax,ebax,ecenterx,ecentery,ephi = fitEllipseContour(ex,ey)
                
                ## Convert back to lat/lon
                elons = np.where(exx<0,np.where(eyy>0,np.arctan(eyy/exx)+math.pi,np.arctan(eyy/exx)-math.pi),np.arctan(eyy/exx))
                yysinxxlon = eyy/np.sin(elons)
                elats = -2*(np.arctan(yysinxxlon) - (math.pi/4.0))
                elats = np.rad2deg(elats)
                elons = np.rad2deg(elons)
                
                ## Still not really sure what this is for...
                for g in range(1,len(elats)):
                    if abs(elats[g]-elats[g-1]) > 1.5:
                        elats[g] = elats[g-1]
                
                ## Center points back to lat/lon
                cenlon = np.where(ecenterx<0,np.where(ecentery>0,np.arctan(ecentery/ecenterx)+math.pi,np.arctan(ecentery/ecenterx)-math.pi),np.arctan(ecentery/ecenterx))
                ysinlon = ecentery/np.sin(cenlon)
                cenlat = np.rad2deg(-2 * (np.arctan(ysinlon) - (math.pi/4.0)))
                cenlon = np.rad2deg(cenlon)
                print("Center of ellipse:",cenlat,"N",cenlon,"E")
                
                ## Calculate endpoints of the axes of the vortex, convert to lat/lon
                xa = eaax * np.cos(ephi)
                ya = eaax * np.sin(ephi)
                xb = ebax * np.sin(ephi)
                yb = ebax * np.cos(ephi)
                endx = np.array([ecenterx+xa,ecenterx-xa,ecenterx+xb,ecenterx-xb])
                endy = np.array([ecentery+ya,ecentery-ya,ecentery-yb,ecentery+yb])
                endlon = np.where(endx < 0,np.where(endy>0,np.arctan(endy/endx)+math.pi,np.arctan(endy/endx)-math.pi),np.arctan(endy/endx))
                ysinelon = endy/np.sin(endlon)
                endlat = np.rad2deg(-2 *(np.arctan(ysinelon) - (math.pi/4)))
                endlon = np.rad2deg(endlon)
                
                ## Calc great circle distances; (still necessary?) 
                a1gc = great_circle((endlat[0],endlon[0]),(cenlat,cenlon)).km
                a2gc = great_circle((cenlat,cenlon),(cenlat,endlat[1])).km
                b1gc = great_circle((endlat[2],endlon[2]),(cenlat,cenlon)).km
                b2gc = great_circle((cenlat,cenlon),(endlat[3],endlon[3])).km
                
                ephi = np.rad2deg(ephi)
                if a1gc < b1gc:
                    ephi -= 90
                    print("Emetrics phi:",ephi)
                if ephi < -45:
                    ephi += 180
                ratio = a1gc/b1gc
                if ratio < 1.0:
                    ratio = 1.0/ratio
                size = math.pi*a1gc*b1gc
                estr = times[t].strftime("%H %Y%m%d%H")+" "+f'{cenlon:.3f}'+" "+f'{cenlat:.3f}'+" "+f'{a1gc:.3f}'+" "+f'{b1gc:.3f}'+" "+f'{ephi:.3f}'+" "+f'{u.values:.3f}'+" "+f'{temp.values:.3f}'
                
                exy=np.array(list(zip(elons,elats)))
                plt.plot(elons,elats,color='red',transform=ccrs.Geodetic(),linewidth=mem_lw)
                plt.plot(cenlon,cenlat,markersize=4,marker='o',color='red',transform=ccrs.PlateCarree())
                #ax.add_patch(Polygon(exy,closed=True,color='magenta',fill=False,lw=mem_lw,transform=ccrs.PlateCarree()))

                ##append values to the empty lists for this day/iteration
                ratio1.append(ratio)
                uvalues1.append(u.values)
                cenlat1.append(cenlat)
                cenlon1.append(cenlon)
                size1.append(size)
                ephi1.append(ephi)

            #check how many vortices and average is >= 1 and inout nan if 0
            if len(ratio1) >= 1:
                rat.append(np.mean(ratio1))
                wind.append(np.mean(uvalues1))
                cenlt.append(np.median(cenlat1))
                cenln.append(np.median(cenlon1))
                sz.append(np.mean(size1))
                ep.append(np.mean(ephi1))
                num.append(len(ratio1))
            
            if len(ratio1) == 0:
                rat.append(np.nan)
                wind.append(np.nan)
                cenlt.append(np.nan)
                cenln.append(np.nan)
                sz.append(np.nan)
                ep.append(np.nan)
                num.append(0)
                
            print("right after savefig")
            plot_label = times[t].strftime("%Y%m%d%H")
            #plt.savefig("./10hPa/"+str(year)+"/10hPa_"+str(plot_label)+".png")
    		#plt.savefig(plots_loc+"ellipse"+str(plot_lev)+"_"+str(fhr/6)+".png",format='png')
            print("right after savefig")
            #image_no += 1
            #plt.close()
            plt.clf()
            
        if plot_lev == 10:
            print('lev 10')
            rat10[:] = rat[:]
            wind10[:] = wind[:]
            cenlat10[:] = cenlt[:]
            cenlon10[:] = cenln[:]
            size10[:] = sz[:]
            ep10[:] = ep[:]
            numb10[:] = num[:]
            
        if plot_lev == 50:
            print('lev 50')
            rat50[:] = rat[:]
            wind50[:] = wind[:]
            cenlat50[:] = cenlt[:]
            cenlon50[:] = cenln[:]
            size50[:] = sz[:]
            ep50[:] = ep[:]
            numb50[:] = num[:]
    
    return rat10,rat50,wind10,wind50,cenlat10,cenlat50,cenlon10,cenlon50,size10,size50,ep10,ep50,numb10,numb50; 

########################################################################################################################
########################################################################################################################

##definition statement for calculating ellipse diagnostics with values from the second ellipse recorded...

def ERA5_two_ellipse_calc(lev_list,contour_list,windlat,temp_lev,y,m,t_d,
                     rat10,wind10,cenlat10,cenlon10,size10,ep10,numb10,
                     rat50,wind50,cenlat50,cenlon50,size50,ep50,numb50):
    #Start with defining desired period of record
    year = y
    month = m
    day = 1
    date1 = dt.datetime(year,month,day,0)  #first date to plot
    total_days = t_d
    first_fhr = 0
    hours = total_days * 24
    hr_inc = 24
    times = [date1 + dt.timedelta(hours=x) for x in range(0,hours,hr_inc)]
    date_list = netCDF4.date2num(times,units="hours since 1800-01-01 00:00:00",calendar="gregorian") #change dates to netcdf times

    ##start loop for calculation and plotting
    for plot_lev,the_contour in zip(lev_list,contour_list): #loop through each level
        print("Opening ERA5 data "+str(plot_lev)+" and contour level "+str(the_contour)  )

        #empty lists for appending diagnostics
        rat = []
        cenlt = []
        cenln = []
        wind = []
        sz = []
        ep = []

        rat2 = []
        cenlt2 = []
        cenln2 = []
        wind2 = []
        sz2 = []
        ep2 = []
        
        num = []
        
        #text_file2 = "./"+date1.strftime("%Y%m%d_%H%M")+"_temp"+str(plot_lev)+".txt"
        #fout2 = open(text_file2,'w')
        if plot_lev == 10:
            ###THESE LINES MUST BE CHANGED FOR WHERE YOUR DATA IS STORED
            gfile = xr.open_dataset("../era5/gph/era5_gph_"+str(year)+".nc")
            g_files = gfile["z"]
            tfile = xr.open_dataset("../era5/t/era5_t_"+str(year)+".nc")
            t_files = tfile["t"]
            ufile = xr.open_dataset("../era5/u/era5_uv_"+str(year)+".nc")
            u_files = ufile["u"]
            
            print("Partitioning and averaging data.")
            print("First, g data")
            g_data = g_files.loc[dict(latitude=slice(90,0),level=lev_list)]
            print("u data")
            u_data = u_files.loc[dict(latitude=slice(windlat+1,windlat-1),level=lev_list)].mean(dim='longitude')
            print("t data")
            t_data = t_files.loc[dict(latitude=slice(90,75),level=lev_list)]
            tlats = t_data['latitude'].values
            weights=np.cos(np.deg2rad(t_data.latitude))
            weights.name = "weights"
            
            print("lat and lon")
            lats_era = g_data["latitude"].values
            lons_era = g_data["longitude"].values
        
        print(f"level: {plot_lev}")
        for date in date_list: #loop through each day 
            print(f"date: {date_list}")
            ind = np.where(date_list == date)[0]
            for d in ind:
                t = d
                height = g_data.loc[dict(time=times[t],level=plot_lev)] ##convert to hPa, gph
                height = height/9.81
                u = u_data.loc[dict(time=times[t],level=plot_lev)].mean(dim='latitude')
                tdata = t_data.loc[dict(time=times[t],level=plot_lev)].mean(dim='longitude')
            temp = tdata.weighted(weights).mean(dim='latitude')
            
            formatted_date = times[t].strftime("%Y%m%d")
            emark = []
            eline = []
            cs_temp = []
            print("##########################")
            print(date)
            print("Date "+times[t].strftime("%H UTC %d %b %Y"))
            valid_label = (times[t].strftime("Valid: %H UTC %d %b %Y"))
            
            plt.Figure(figsize=(15,15),dpi=120)#figure(figsize=(12,12),dpi=1200) <---Set fig size here!
            #ax = plt.subplot(1,1,1,projection=ccrs.Orthographic(0,90))
            ax = plt.axes(projection=ccrs.Orthographic(0,90))
            ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='dimgray',facecolor='none')
            #ax.set_extent([-180,180,0,90],ccrs.PlateCarree())
            ax.outline_patch.set_edgecolor('none')
            gl = ax.gridlines(color="grey",linestyle=":",linewidth=0.5)#=lat_lon_width)
            
            clevs = range(18000,33500,250)
            
            [x,y] = np.meshgrid(lons_era,lats_era)
            
            plt.title(str(plot_lev)+"hPa Elliptical Diagnostics\n"+valid_label)
            
            mem_color = "blue"
            mem_lw = 3.0
            mem_ms = 8
            cont_color = "#aec2e2"
            cont_lw = 1.0
            
            cs_temp.append(ax.contour(x,y,height,levels=clevs,linewidths=cont_lw,colors=cont_color))
            
            line_w = np.array(clevs,dtype=np.float)
            for i,c in enumerate(clevs):
                if c == the_contour:
                    line_w[i] = 1.5 # Make the contour we're using for the ellipse calc fatter
                else:
                    line_w[i] = 0.5
            print(clevs)
            ax.contour(x,y,height,clevs,transform=ccrs.PlateCarree(),extend='both',colors='black',linewidths=line_w)
            #check number/presence of contours
            try:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                print("in try")
                print("Number of ellipses: ",len(isoline_list))
            except:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []
            
            if len(isoline_list) > 0:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                #print(isoline_list)
                print("in ellipse ...")
                #quit()
                print("Number of ellipses: ",len(isoline_list))
            else:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []

            ##number of contours
            isocount = 0
            small = 0
            #empty lists for upcoming metric calculations
            ratio1 = []
            uvalues1 = []
            cenlat1 = []
            cenlon1 = []
            size1 = []
            ephi1 = []
            
            for isoline in isoline_list:
                ##next few lines will check if the contour crossing the prime meridian is one continuous contour or two separate ones
                #print('isoline',isoline)
                #[iso_lon,iso_lat] = mm(isoline[:,0],isoline[:,1],inverse=True)
                [iso_lon,iso_lat] = [isoline[:,0],isoline[:,1]]
                if len(iso_lon)<15:
                    print("-----Not analyzing ellipse with",len(iso_lon),"points, continuing...")  # Check for size!
                    small = 1
                    continue
                # Checking to see if contours are closed (0) or need to be joined (1)  Checking before convert lat/lon to radians
                lon_diff = abs(iso_lon[0] - iso_lon[len(iso_lon)-1])
                lat_diff = abs(iso_lat[0] - iso_lat[len(iso_lat)-1])
                join = 0
                if lon_diff > 1 or lat_diff > 1:
                    join = 1
                    print("Diffs lat/lon: ",lat_diff,lon_diff)
                iso_lon = np.deg2rad(iso_lon)
                iso_lat = np.deg2rad(iso_lat)
                
                ex = np.array((np.cos(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ey = np.array((np.sin(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ## does this contour include the pole?
                overpole = point_inside_polygon(0,0,ex,ey)  #returns true if poly includes the pole, false if not 
                print("It is",overpole, "that the contour includes the pole")
                if not overpole and len(isoline_list) > 1 and join > 0 and small < 1:
                    if isocount > 0:
                        ex = np.append(ex,ex2)
                        ey = np.append(ey,ey2)
                        print("Point didn't include pole - add to these to previous set")
                    else:
                        ex2 = ex
                        ey2 = ey
                        print("Point didn't include pole - keeping these to add to next set")
                        isocount = 1
                        continue
                ##begin actual calculation        
                print("Running ellipse diagnostic now")
                exx,eyy,eaax,ebax,ecenterx,ecentery,ephi = fitEllipseContour(ex,ey)
                
                ## Convert back to lat/lon
                elons = np.where(exx<0,np.where(eyy>0,np.arctan(eyy/exx)+math.pi,np.arctan(eyy/exx)-math.pi),np.arctan(eyy/exx))
                yysinxxlon = eyy/np.sin(elons)
                elats = -2*(np.arctan(yysinxxlon) - (math.pi/4.0))
                elats = np.rad2deg(elats)
                elons = np.rad2deg(elons)
                
                ## Still not really sure what this is for...
                for g in range(1,len(elats)):
                    if abs(elats[g]-elats[g-1]) > 1.5:
                        elats[g] = elats[g-1]
                
                ## Center points back to lat/lon
                cenlon = np.where(ecenterx<0,np.where(ecentery>0,np.arctan(ecentery/ecenterx)+math.pi,np.arctan(ecentery/ecenterx)-math.pi),np.arctan(ecentery/ecenterx))
                ysinlon = ecentery/np.sin(cenlon)
                cenlat = np.rad2deg(-2 * (np.arctan(ysinlon) - (math.pi/4.0)))
                cenlon = np.rad2deg(cenlon)
                print("Center of ellipse:",cenlat,"N",cenlon,"E")
                
                ## Calculate endpoints of the axes of the vortex, convert to lat/lon
                xa = eaax * np.cos(ephi)
                ya = eaax * np.sin(ephi)
                xb = ebax * np.sin(ephi)
                yb = ebax * np.cos(ephi)
                endx = np.array([ecenterx+xa,ecenterx-xa,ecenterx+xb,ecenterx-xb])
                endy = np.array([ecentery+ya,ecentery-ya,ecentery-yb,ecentery+yb])
                endlon = np.where(endx < 0,np.where(endy>0,np.arctan(endy/endx)+math.pi,np.arctan(endy/endx)-math.pi),np.arctan(endy/endx))
                ysinelon = endy/np.sin(endlon)
                endlat = np.rad2deg(-2 *(np.arctan(ysinelon) - (math.pi/4)))
                endlon = np.rad2deg(endlon)
                
                ## Calc great circle distances; (still necessary?) 
                a1gc = great_circle((endlat[0],endlon[0]),(cenlat,cenlon)).km
                a2gc = great_circle((cenlat,cenlon),(cenlat,endlat[1])).km
                b1gc = great_circle((endlat[2],endlon[2]),(cenlat,cenlon)).km
                b2gc = great_circle((cenlat,cenlon),(endlat[3],endlon[3])).km
                
                ephi = np.rad2deg(ephi)
                if a1gc < b1gc:
                    ephi -= 90
                    print("Emetrics phi:",ephi)
                if ephi < -45:
                    ephi += 180
                ratio = a1gc/b1gc
                if ratio < 1.0:
                    ratio = 1.0/ratio
                size = math.pi*a1gc*b1gc
                estr = times[t].strftime("%H %Y%m%d%H")+" "+f'{cenlon:.3f}'+" "+f'{cenlat:.3f}'+" "+f'{a1gc:.3f}'+" "+f'{b1gc:.3f}'+" "+f'{ephi:.3f}'+" "+f'{u.values:.3f}'+" "+f'{temp.values:.3f}'
                
                exy=np.array(list(zip(elons,elats)))
                plt.plot(elons,elats,color='red',transform=ccrs.Geodetic(),linewidth=mem_lw)
                plt.plot(cenlon,cenlat,markersize=4,marker='o',color='red',transform=ccrs.PlateCarree())
                #ax.add_patch(Polygon(exy,closed=True,color='magenta',fill=False,lw=mem_lw,transform=ccrs.PlateCarree()))

                ##append values to the empty lists for this day/iteration
                ratio1.append(ratio)
                uvalues1.append(u.values)
                cenlat1.append(cenlat)
                cenlon1.append(cenlon)
                size1.append(size)
                ephi1.append(ephi)

            #check how many vortices and average is >= 1 and inout nan if 0
            if len(ratio1) > 1:
                rat.append(ratio1[0])
                wind.append(uvalues1[0])
                cenlt.append(cenlat1[0])
                cenln.append(cenlon1[0])
                sz.append(size1[0])
                ep.append(ephi1[0])

                rat2.append(ratio1[1])
                wind2.append(uvalues1[1])
                cenlt2.append(cenlat1[1])
                cenln2.append(cenlon1[1])
                sz2.append(size1[1])
                ep2.append(ephi1[1])
                
                num.append(len(ratio1))
                
            if len(ratio1) == 1:
                rat.append(ratio1[0])
                wind.append(uvalues1[0])
                cenlt.append(cenlat1[0])
                cenln.append(cenlon1[0])
                sz.append(size1[0])
                ep.append(ephi1[0])

                rat2.append(np.nan)
                wind2.append(np.nan)
                cenlt2.append(np.nan)
                cenln2.append(np.nan)
                sz2.append(np.nan)
                ep2.append(np.nan)
                
                num.append(len(ratio1))
            
            if len(ratio1) == 0:
                rat.append(np.nan)
                wind.append(np.nan)
                cenlt.append(np.nan)
                cenln.append(np.nan)
                sz.append(np.nan)
                ep.append(np.nan)

                rat2.append(np.nan)
                wind2.append(np.nan)
                cenlt2.append(np.nan)
                cenln2.append(np.nan)
                sz2.append(np.nan)
                ep2.append(np.nan)
                
                num.append(0)
                
            print("right after savefig")
            plot_label = times[t].strftime("%Y%m%d%H")
            #plt.savefig("./10hPa/"+str(year)+"/10hPa_"+str(plot_label)+".png")
    		#plt.savefig(plots_loc+"ellipse"+str(plot_lev)+"_"+str(fhr/6)+".png",format='png')
            print("right after savefig")
            #image_no += 1
            #plt.close()
            plt.clf()
            
        if plot_lev == 10:
            print('lev 10')
            rat10[:] = [rat[:],rat2[:]]
            wind10[:] = [wind[:],wind2[:]]
            cenlat10[:] = [cenlt[:],cenlt2[:]]
            cenlon10[:] = [cenln[:],cenln2[:]]
            size10[:] = [sz[:],sz2[:]]
            ep10[:] = [ep[:],ep2[:]]
            numb10[:] = num[:]
            
        if plot_lev == 50:
            print('lev 50')
            rat50[:] = [rat[:],rat2[:]]
            wind50[:] = [wind[:],wind2[:]]
            cenlat50[:] = [cenlt[:],cenlt2[:]]
            cenlon50[:] = [cenln[:],cenln2[:]]
            size50[:] = [sz[:],sz2[:]]
            ep50[:] = [ep[:],ep2[:]]
            numb50[:] = num[:]
    
    return rat10,rat50,wind10,wind50,cenlat10,cenlat50,cenlon10,cenlon50,size10,size50,ep10,ep50,numb10,numb50; 

########################################################################################################################
########################################################################################################################

def ERAI_one_ellipse_calc(lev_list,contour_list,windlat,temp_lev,y,m,t_d,
                     rat10,wind10,cenlat10,cenlon10,size10,ep10,numb10,
                     rat50,wind50,cenlat50,cenlon50,size50,ep50,numb50):
       #Start with defining desired period of record
    year = y
    month = m
    day = 1
    date1 = dt.datetime(year,month,day,0)  #first date to plot
    total_days = t_d
    first_fhr = 0
    hours = total_days * 24
    hr_inc = 24
    times = [date1 + dt.timedelta(hours=x) for x in range(0,hours,hr_inc)]
    date_list = netCDF4.date2num(times,units="hours since 1800-01-01 00:00:00",calendar="gregorian") #change dates to netcdf times

    ##start loop for calculation and plotting
    for plot_lev,the_contour in zip(lev_list,contour_list): #loop through each level
        print("Opening ERAI data "+str(plot_lev)+" and contour level "+str(the_contour)  )

        #empty lists for appending diagnostics
        rat = []
        cenlt = []
        cenln = []
        wind = []
        sz = []
        ep = []
        
        num = []
        
        #text_file2 = "./"+date1.strftime("%Y%m%d_%H%M")+"_temp"+str(plot_lev)+".txt"
        #fout2 = open(text_file2,'w')
        if plot_lev == 10:
            ###THESE LINES MUST BE CHANGED FOR WHERE YOUR DATA IS STORED
            gfile = xr.open_dataset("/erai/"+str(year)+"/g."+str(year)+".nc")
            g_files = gfile["g"]
            tfile = xr.open_dataset("/erai/"+str(year)+"/t."+str(year)+".nc")
            t_files = tfile["t"]
            ufile = xr.open_dataset("/erai/"+str(year)+"/u."+str(year)+".nc")
            u_files = ufile["u"]
            
            g_data = g_files.loc[dict(lat=slice(0,90),lev=lev_list)]
            u_data = u_files.loc[dict(lat=slice(windlat-1,windlat+1),lev=lev_list)].mean(dim='lon')
            t_data = t_files.loc[dict(lat=slice(75,90),lev=lev_list)]
            tlats = t_data['lat'].values
            weights=np.cos(np.deg2rad(t_data.lat))
            weights.name = "weights"
            
            lats_era = g_data["lat"].values
            lons_era = g_data["lon"].values
        
        for date in date_list: #loop through each day 
            #print(f"date: {date_list}")
            ind = np.where(date_list == date)[0]
            for d in ind:
                t = d
                height = g_data.loc[dict(time=times[t],lev=plot_lev)]
                u = u_data.loc[dict(time=times[t],lev=plot_lev)].mean(dim='lat')
                tdata = t_data.loc[dict(time=times[t],lev=plot_lev)].mean(dim='lon')
            temp = tdata.weighted(weights).mean(dim='lat')
            
            formatted_date = times[t].strftime("%Y%m%d")
            emark = []
            eline = []
            cs_temp = []
            print("##########################")
            print(date)
            print("Date "+times[t].strftime("%H UTC %d %b %Y"))
            valid_label = (times[t].strftime("Valid: %H UTC %d %b %Y"))
            
            plt.Figure(figsize=(15,15),dpi=120)#figure(figsize=(12,12),dpi=1200) <---Set fig size here!
            #ax = plt.subplot(1,1,1,projection=ccrs.Orthographic(0,90))
            ax = plt.axes(projection=ccrs.Orthographic(0,90))
            ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='dimgray',facecolor='none')
            #ax.set_extent([-180,180,0,90],ccrs.PlateCarree())
            ax.outline_patch.set_edgecolor('none')
            gl = ax.gridlines(color="grey",linestyle=":",linewidth=0.5)#=lat_lon_width)
            
            clevs = range(18000,33500,250)
            
            [x,y] = np.meshgrid(lons_era,lats_era)
            
            plt.title(str(plot_lev)+"hPa Elliptical Diagnostics\n"+valid_label)
            
            mem_color = "blue"
            mem_lw = 3.0
            mem_ms = 8
            cont_color = "#aec2e2"
            cont_lw = 1.0
            
            cs_temp.append(ax.contour(x,y,height,levels=clevs,linewidths=cont_lw,colors=cont_color))
            
            line_w = np.array(clevs,dtype=np.float)
            for i,c in enumerate(clevs):
                if c == the_contour:
                    line_w[i] = 1.5 # Make the contour we're using for the ellipse calc fatter
                else:
                    line_w[i] = 0.5
            print(clevs)
            ax.contour(x,y,height,clevs,transform=ccrs.PlateCarree(),extend='both',colors='black',linewidths=line_w)
            #check number/presence of contours
            try:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                print("in try")
                print("Number of ellipses: ",len(isoline_list))
            except:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []
            
            if len(isoline_list) > 0:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                #print(isoline_list)
                print("in ellipse ...")
                #quit()
                print("Number of ellipses: ",len(isoline_list))
            else:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []

            ##number of contours
            isocount = 0
            small = 0
            #empty lists for upcoming metric calculations
            ratio1 = []
            uvalues1 = []
            cenlat1 = []
            cenlon1 = []
            size1 = []
            ephi1 = []
            
            for isoline in isoline_list:
                ##next few lines will check if the contour crossing the prime meridian is one continuous contour or two separate ones
                #print('isoline',isoline)
                #[iso_lon,iso_lat] = mm(isoline[:,0],isoline[:,1],inverse=True)
                [iso_lon,iso_lat] = [isoline[:,0],isoline[:,1]]
                if len(iso_lon)<15:
                    print("-----Not analyzing ellipse with",len(iso_lon),"points, continuing...")  # Check for size!
                    small = 1
                    continue
                # Checking to see if contours are closed (0) or need to be joined (1)  Checking before convert lat/lon to radians
                lon_diff = abs(iso_lon[0] - iso_lon[len(iso_lon)-1])
                lat_diff = abs(iso_lat[0] - iso_lat[len(iso_lat)-1])
                join = 0
                if lon_diff > 1 or lat_diff > 1:
                    join = 1
                    print("Diffs lat/lon: ",lat_diff,lon_diff)
                iso_lon = np.deg2rad(iso_lon)
                iso_lat = np.deg2rad(iso_lat)
                
                ex = np.array((np.cos(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ey = np.array((np.sin(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ## does this contour include the pole?
                overpole = point_inside_polygon(0,0,ex,ey)  #returns true if poly includes the pole, false if not 
                print("It is",overpole, "that the contour includes the pole")
                if not overpole and len(isoline_list) > 1 and join > 0 and small < 1:
                    if isocount > 0:
                        ex = np.append(ex,ex2)
                        ey = np.append(ey,ey2)
                        print("Point didn't include pole - add to these to previous set")
                    else:
                        ex2 = ex
                        ey2 = ey
                        print("Point didn't include pole - keeping these to add to next set")
                        isocount = 1
                        continue
                ##begin actual calculation        
                print("Running ellipse diagnostic now")
                exx,eyy,eaax,ebax,ecenterx,ecentery,ephi = fitEllipseContour(ex,ey)
                
                ## Convert back to lat/lon
                elons = np.where(exx<0,np.where(eyy>0,np.arctan(eyy/exx)+math.pi,np.arctan(eyy/exx)-math.pi),np.arctan(eyy/exx))
                yysinxxlon = eyy/np.sin(elons)
                elats = -2*(np.arctan(yysinxxlon) - (math.pi/4.0))
                elats = np.rad2deg(elats)
                elons = np.rad2deg(elons)
                
                ## Still not really sure what this is for...
                for g in range(1,len(elats)):
                    if abs(elats[g]-elats[g-1]) > 1.5:
                        elats[g] = elats[g-1]
                
                ## Center points back to lat/lon
                cenlon = np.where(ecenterx<0,np.where(ecentery>0,np.arctan(ecentery/ecenterx)+math.pi,np.arctan(ecentery/ecenterx)-math.pi),np.arctan(ecentery/ecenterx))
                ysinlon = ecentery/np.sin(cenlon)
                cenlat = np.rad2deg(-2 * (np.arctan(ysinlon) - (math.pi/4.0)))
                cenlon = np.rad2deg(cenlon)
                print("Center of ellipse:",cenlat,"N",cenlon,"E")
                
                ## Calculate endpoints of the axes of the vortex, convert to lat/lon
                xa = eaax * np.cos(ephi)
                ya = eaax * np.sin(ephi)
                xb = ebax * np.sin(ephi)
                yb = ebax * np.cos(ephi)
                endx = np.array([ecenterx+xa,ecenterx-xa,ecenterx+xb,ecenterx-xb])
                endy = np.array([ecentery+ya,ecentery-ya,ecentery-yb,ecentery+yb])
                endlon = np.where(endx < 0,np.where(endy>0,np.arctan(endy/endx)+math.pi,np.arctan(endy/endx)-math.pi),np.arctan(endy/endx))
                ysinelon = endy/np.sin(endlon)
                endlat = np.rad2deg(-2 *(np.arctan(ysinelon) - (math.pi/4)))
                endlon = np.rad2deg(endlon)
                
                ## Calc great circle distances; (still necessary?) 
                a1gc = great_circle((endlat[0],endlon[0]),(cenlat,cenlon)).km
                a2gc = great_circle((cenlat,cenlon),(cenlat,endlat[1])).km
                b1gc = great_circle((endlat[2],endlon[2]),(cenlat,cenlon)).km
                b2gc = great_circle((cenlat,cenlon),(endlat[3],endlon[3])).km
                
                ephi = np.rad2deg(ephi)
                if a1gc < b1gc:
                    ephi -= 90
                    print("Emetrics phi:",ephi)
                if ephi < -45:
                    ephi += 180
                ratio = a1gc/b1gc
                if ratio < 1.0:
                    ratio = 1.0/ratio
                size = math.pi*a1gc*b1gc
                estr = times[t].strftime("%H %Y%m%d%H")+" "+f'{cenlon:.3f}'+" "+f'{cenlat:.3f}'+" "+f'{a1gc:.3f}'+" "+f'{b1gc:.3f}'+" "+f'{ephi:.3f}'+" "+f'{u.values:.3f}'+" "+f'{temp.values:.3f}'
                
                exy=np.array(list(zip(elons,elats)))
                plt.plot(elons,elats,color='red',transform=ccrs.Geodetic(),linewidth=mem_lw)
                plt.plot(cenlon,cenlat,markersize=4,marker='o',color='red',transform=ccrs.PlateCarree())
                #ax.add_patch(Polygon(exy,closed=True,color='magenta',fill=False,lw=mem_lw,transform=ccrs.PlateCarree()))

                ##append values to the empty lists for this day/iteration
                ratio1.append(ratio)
                uvalues1.append(u.values)
                cenlat1.append(cenlat)
                cenlon1.append(cenlon)
                size1.append(size)
                ephi1.append(ephi)
                
            #check how many vortices and average is >= 1 and inout nan if 0
            if len(ratio1) >= 1:
                rat.append(np.mean(ratio1))
                wind.append(np.mean(uvalues1))
                cenlt.append(np.median(cenlat1))
                cenln.append(np.median(cenlon1))
                sz.append(np.mean(size1))
                ep.append(np.mean(ephi1))
                num.append(len(ratio1))
            
            if len(ratio1) == 0:
                rat.append(np.nan)
                wind.append(np.nan)
                cenlt.append(np.nan)
                cenln.append(np.nan)
                sz.append(np.nan)
                ep.append(np.nan)
                num.append(0)
                
            print("right after savefig")
            plot_label = times[t].strftime("%Y%m%d%H")
            #plt.savefig("./10hPa/"+str(year)+"/10hPa_"+str(plot_label)+".png")
    		#plt.savefig(plots_loc+"ellipse"+str(plot_lev)+"_"+str(fhr/6)+".png",format='png')
            print("right after savefig")
            #image_no += 1
            #plt.close()
            plt.clf()
            
        if plot_lev == 10:
            print('lev 10')
            rat10[:] = rat[:]
            wind10[:] = wind[:]
            cenlat10[:] = cenlt[:]
            cenlon10[:] = cenln[:]
            size10[:] = sz[:]
            ep10[:] = ep[:]
            numb10[:] = num[:]
            
        if plot_lev == 50:
            print('lev 50')
            rat50[:] = rat[:]
            wind50[:] = wind[:]
            cenlat50[:] = cenlt[:]
            cenlon50[:] = cenln[:]
            size50[:] =sz[:]
            ep50[:] = ep[:]
            numb50[:] = num[:]
    
    return rat10,rat50,wind10,wind50,cenlat10,cenlat50,cenlon10,cenlon50,size10,size50,ep10,ep50,numb10,numb50; 

########################################################################################################################
########################################################################################################################
def ERAI_two_ellipse_calc(lev_list,contour_list,windlat,temp_lev,y,m,t_d,
                     rat10,wind10,cenlat10,cenlon10,size10,ep10,numb10,
                     rat50,wind50,cenlat50,cenlon50,size50,ep50,numb50):
	#Start with defining desired period of record
    year = y
    month = m
    day = 1
    date1 = dt.datetime(year,month,day,0)  #first date to plot
    total_days = t_d
    first_fhr = 0
    hours = total_days * 24
    hr_inc = 24
    times = [date1 + dt.timedelta(hours=x) for x in range(0,hours,hr_inc)]
    date_list = netCDF4.date2num(times,units="hours since 1800-01-01 00:00:00",calendar="gregorian") #change dates to netcdf times
    
    for plot_lev,the_contour in zip(lev_list,contour_list): #loop through each level
        print("Opening ERAI data "+str(plot_lev)+" and contour level "+str(the_contour))
        
        rat = []
        cenlt = []
        cenln = []
        wind = []
        sz = []
        ep = []

        rat2 = []
        cenlt2 = []
        cenln2 = []
        wind2 = []
        sz2 = []
        ep2 = []
        
        num = []
        
        if plot_lev == 10:
            ###THESE LINES MUST BE CHANGED FOR WHERE YOUR DATA IS STORED
            gfile = xr.open_dataset("/erai/"+str(year)+"/g."+str(year)+".nc")
            g_files = gfile["g"]
            tfile = xr.open_dataset("/erai/"+str(year)+"/t."+str(year)+".nc")
            t_files = tfile["t"]
            ufile = xr.open_dataset("/erai/"+str(year)+"/u."+str(year)+".nc")
            u_files = ufile["u"]
            
            g_data = g_files.loc[dict(lat=slice(0,90),lev=lev_list)]
            u_data = u_files.loc[dict(lat=slice(windlat-1,windlat+1),lev=lev_list)].mean(dim='lon')
            t_data = t_files.loc[dict(lat=slice(75,90),lev=lev_list)]
            tlats = t_data['lat'].values
            weights=np.cos(np.deg2rad(t_data.lat))
            weights.name = "weights"
            
            lats_era = g_data["lat"].values
            lons_era = g_data["lon"].values
        
        print(f"level: {plot_lev}")
        for date in date_list: #loop through each day 
            #print(f"date: {date_list}")
            ind = np.where(date_list == date)[0]
            for d in ind:
                t = d
                height = g_data.loc[dict(time=times[t],lev=plot_lev)]
                u = u_data.loc[dict(time=times[t],lev=plot_lev)].mean(dim='lat')
                tdata = t_data.loc[dict(time=times[t],lev=plot_lev)].mean(dim='lon')
            temp = tdata.weighted(weights).mean(dim='lat')
            
            formatted_date = times[t].strftime("%Y%m%d")
            emark = []
            eline = []
            cs_temp = []
            print("##########################")
            print(date)
            print("Date "+times[t].strftime("%H UTC %d %b %Y"))
            valid_label = (times[t].strftime("Valid: %H UTC %d %b %Y"))
            
            plt.Figure(figsize=(15,15),dpi=120)#figure(figsize=(12,12),dpi=1200) <---Set fig size here!
            #ax = plt.subplot(1,1,1,projection=ccrs.Orthographic(0,90))
            ax = plt.axes(projection=ccrs.Orthographic(0,90))
            ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='dimgray',facecolor='none')
            #ax.set_extent([-180,180,0,90],ccrs.PlateCarree())
            ax.outline_patch.set_edgecolor('none')
            gl = ax.gridlines(color="grey",linestyle=":",linewidth=0.5)#=lat_lon_width)
            
            clevs = range(18000,33500,250)
            
            [x,y] = np.meshgrid(lons_era,lats_era)
            
            plt.title(str(plot_lev)+"hPa Elliptical Diagnostics\n"+valid_label)
            
            mem_color = "blue"
            mem_lw = 3.0
            mem_ms = 8
            cont_color = "#aec2e2"
            cont_lw = 1.0
            
            cs_temp.append(ax.contour(x,y,height,levels=clevs,linewidths=cont_lw,colors=cont_color))
            
            line_w = np.array(clevs,dtype=np.float)
            for i,c in enumerate(clevs):
                if c == the_contour:
                    line_w[i] = 1.5 # Make the contour we're using for the ellipse calc fatter
                else:
                    line_w[i] = 0.5
            print(clevs)
            ax.contour(x,y,height,clevs,transform=ccrs.PlateCarree(),extend='both',colors='black',linewidths=line_w)
            #check number/presence of contours
            try:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                print("in try")
                print("Number of ellipses: ",len(isoline_list))
            except:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []
            
            if len(isoline_list) > 0:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                #print(isoline_list)
                print("in ellipse ...")
                #quit()
                print("Number of ellipses: ",len(isoline_list))
            else:
                print("Contours not found for",the_contour,"meter level at",plot_lev,"hPa.")
                isoline_list = []

            ##number of contours
            isocount = 0
            small = 0
            #empty lists for upcoming metric calculations
            ratio1 = []
            uvalues1 = []
            cenlat1 = []
            cenlon1 = []
            size1 = []
            ephi1 = []
            
            for isoline in isoline_list:
                ##next few lines will check if the contour crossing the prime meridian is one continuous contour or two separate ones
                #print('isoline',isoline)
                #[iso_lon,iso_lat] = mm(isoline[:,0],isoline[:,1],inverse=True)
                [iso_lon,iso_lat] = [isoline[:,0],isoline[:,1]]
                if len(iso_lon)<15:
                    print("-----Not analyzing ellipse with",len(iso_lon),"points, continuing...")  # Check for size!
                    small = 1
                    continue
                # Checking to see if contours are closed (0) or need to be joined (1)  Checking before convert lat/lon to radians
                lon_diff = abs(iso_lon[0] - iso_lon[len(iso_lon)-1])
                lat_diff = abs(iso_lat[0] - iso_lat[len(iso_lat)-1])
                join = 0
                if lon_diff > 1 or lat_diff > 1:
                    join = 1
                    print("Diffs lat/lon: ",lat_diff,lon_diff)
                iso_lon = np.deg2rad(iso_lon)
                iso_lat = np.deg2rad(iso_lat)
                
                ex = np.array((np.cos(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ey = np.array((np.sin(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ## does this contour include the pole?
                overpole = point_inside_polygon(0,0,ex,ey)  #returns true if poly includes the pole, false if not 
                print("It is",overpole, "that the contour includes the pole")
                if not overpole and len(isoline_list) > 1 and join > 0 and small < 1:
                    if isocount > 0:
                        ex = np.append(ex,ex2)
                        ey = np.append(ey,ey2)
                        print("Point didn't include pole - add to these to previous set")
                    else:
                        ex2 = ex
                        ey2 = ey
                        print("Point didn't include pole - keeping these to add to next set")
                        isocount = 1
                        continue
                ##begin actual calculation        
                print("Running ellipse diagnostic now")
                exx,eyy,eaax,ebax,ecenterx,ecentery,ephi = fitEllipseContour(ex,ey)
                
                ## Convert back to lat/lon
                elons = np.where(exx<0,np.where(eyy>0,np.arctan(eyy/exx)+math.pi,np.arctan(eyy/exx)-math.pi),np.arctan(eyy/exx))
                yysinxxlon = eyy/np.sin(elons)
                elats = -2*(np.arctan(yysinxxlon) - (math.pi/4.0))
                elats = np.rad2deg(elats)
                elons = np.rad2deg(elons)
                
                ## Still not really sure what this is for...
                for g in range(1,len(elats)):
                    if abs(elats[g]-elats[g-1]) > 1.5:
                        elats[g] = elats[g-1]
                
                ## Center points back to lat/lon
                cenlon = np.where(ecenterx<0,np.where(ecentery>0,np.arctan(ecentery/ecenterx)+math.pi,np.arctan(ecentery/ecenterx)-math.pi),np.arctan(ecentery/ecenterx))
                ysinlon = ecentery/np.sin(cenlon)
                cenlat = np.rad2deg(-2 * (np.arctan(ysinlon) - (math.pi/4.0)))
                cenlon = np.rad2deg(cenlon)
                print("Center of ellipse:",cenlat,"N",cenlon,"E")
                
                ## Calculate endpoints of the axes of the vortex, convert to lat/lon
                xa = eaax * np.cos(ephi)
                ya = eaax * np.sin(ephi)
                xb = ebax * np.sin(ephi)
                yb = ebax * np.cos(ephi)
                endx = np.array([ecenterx+xa,ecenterx-xa,ecenterx+xb,ecenterx-xb])
                endy = np.array([ecentery+ya,ecentery-ya,ecentery-yb,ecentery+yb])
                endlon = np.where(endx < 0,np.where(endy>0,np.arctan(endy/endx)+math.pi,np.arctan(endy/endx)-math.pi),np.arctan(endy/endx))
                ysinelon = endy/np.sin(endlon)
                endlat = np.rad2deg(-2 *(np.arctan(ysinelon) - (math.pi/4)))
                endlon = np.rad2deg(endlon)
                
                ## Calc great circle distances; (still necessary?) 
                a1gc = great_circle((endlat[0],endlon[0]),(cenlat,cenlon)).km
                a2gc = great_circle((cenlat,cenlon),(cenlat,endlat[1])).km
                b1gc = great_circle((endlat[2],endlon[2]),(cenlat,cenlon)).km
                b2gc = great_circle((cenlat,cenlon),(endlat[3],endlon[3])).km
                
                ephi = np.rad2deg(ephi)
                if a1gc < b1gc:
                    ephi -= 90
                    print("Emetrics phi:",ephi)
                if ephi < -45:
                    ephi += 180
                ratio = a1gc/b1gc
                if ratio < 1.0:
                    ratio = 1.0/ratio
                size = math.pi*a1gc*b1gc
                estr = times[t].strftime("%H %Y%m%d%H")+" "+f'{cenlon:.3f}'+" "+f'{cenlat:.3f}'+" "+f'{a1gc:.3f}'+" "+f'{b1gc:.3f}'+" "+f'{ephi:.3f}'+" "+f'{u.values:.3f}'+" "+f'{temp.values:.3f}'
                
                exy=np.array(list(zip(elons,elats)))
                plt.plot(elons,elats,color='red',transform=ccrs.Geodetic(),linewidth=mem_lw)
                plt.plot(cenlon,cenlat,markersize=4,marker='o',color='red',transform=ccrs.PlateCarree())
                #ax.add_patch(Polygon(exy,closed=True,color='magenta',fill=False,lw=mem_lw,transform=ccrs.PlateCarree()))

                ##append values to the empty lists for this day/iteration
                ratio1.append(ratio)
                uvalues1.append(u.values)
                cenlat1.append(cenlat)
                cenlon1.append(cenlon)
                size1.append(size)
                ephi1.append(ephi)

            #check how many vortices and average is >= 1 and inout nan if 0
            if len(ratio1) > 1:
                rat.append(ratio1[0])
                wind.append(uvalues1[0])
                cenlt.append(cenlat1[0])
                cenln.append(cenlon1[0])
                sz.append(size1[0])
                ep.append(ephi1[0])

                rat2.append(ratio1[1])
                wind2.append(uvalues1[1])
                cenlt2.append(cenlat1[1])
                cenln2.append(cenlon1[1])
                sz2.append(size1[1])
                ep2.append(ephi1[1])
                
                num.append(len(ratio1))
                
            if len(ratio1) == 1:
                rat.append(ratio1[0])
                wind.append(uvalues1[0])
                cenlt.append(cenlat1[0])
                cenln.append(cenlon1[0])
                sz.append(size1[0])
                ep.append(ephi1[0])

                rat2.append(np.nan)
                wind2.append(np.nan)
                cenlt2.append(np.nan)
                cenln2.append(np.nan)
                sz2.append(np.nan)
                ep2.append(np.nan)
                
                num.append(len(ratio1))
            
            if len(ratio1) == 0:
                rat.append(np.nan)
                wind.append(np.nan)
                cenlt.append(np.nan)
                cenln.append(np.nan)
                sz.append(np.nan)
                ep.append(np.nan)

                rat2.append(np.nan)
                wind2.append(np.nan)
                cenlt2.append(np.nan)
                cenln2.append(np.nan)
                sz2.append(np.nan)
                ep2.append(np.nan)
                
                num.append(0)
                
            print("right after savefig")
            plot_label = times[t].strftime("%Y%m%d%H")
            #plt.savefig("./10hPa/"+str(year)+"/10hPa_"+str(plot_label)+".png")
    		#plt.savefig(plots_loc+"ellipse"+str(plot_lev)+"_"+str(fhr/6)+".png",format='png')
            print("right after savefig")
            #image_no += 1
            #plt.close()
            plt.clf()
            
        if plot_lev == 10:
            print('lev 10')
            rat10[:] = [rat[:],rat2[:]]
            wind10[:] = [wind[:],wind2[:]]
            cenlat10[:] = [cenlt[:],cenlt2[:]]
            cenlon10[:] = [cenln[:],cenln2[:]]
            size10[:] = [sz[:],sz2[:]]
            ep10[:] = [ep[:],ep2[:]]
            numb10[:] = num[:]
            
        if plot_lev == 50:
            print('lev 10')
            rat50[:] = [rat[:],rat2[:]]
            wind50[:] = [wind[:],wind2[:]]
            cenlat50[:] = [cenlt[:],cenlt2[:]]
            cenlon50[:] = [cenln[:],cenln2[:]]
            size50[:] = [sz[:],sz2[:]]
            ep50[:] = [ep[:],ep2[:]]
            numb50[:] = num[:]
    
    return rat10,rat50,wind10,wind50,cenlat10,cenlat50,cenlon10,cenlon50,size10,size50,ep10,ep50,numb10,numb50; 

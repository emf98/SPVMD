import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic
import xarray as xr
import math
from get_ellipse_metrics import get_emetrics_max_min
from fitEllipse2_new import fitEllipseContour
from geopy.distance import great_circle
from matplotlib.patches import Polygon
import copy
import os
os.environ['DISPLAY']='harlan.atmos.albany.edu:1.0'

#### Define functions first:
# determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This is taken from http://www.ariel.com.au/a/python-point-int-poly.html

def point_inside_polygon(x,y,xarr,yarr):

    n = len(xarr)
    inside =False

    p1x = xarr[0]
    p1y = yarr[0]
    for i in range(n+1):
        p2x,p2y = xarr[i % n],yarr[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside;



lev_list = [10,50,30]
contour_list = [30000,20000,23000]

for plot_lev,the_contour in zip(lev_list,contour_list):

    #plot_lev = 10
    ## This must have a time component, or else fhr_inc<24h won't work.
    today = dt.datetime.utcnow()-dt.timedelta(days=0)   ### changed this to 0 day delay
    today = today.replace(hour=0,minute=0,second=0,microsecond=0)
    day_str = today.strftime("%Y%m%d")
    day_label = today.strftime("Init: %a %b %d %Y %H%M UTC")
    first_fhr = 0
    last_fhr = 384
    fhr_inc = 6
    #the_contour = 30000
    text_file = "./plots/emetrics_"+today.strftime("%Y%m%d_%H%M")+"_"+str(plot_lev)
    plots_loc = "./plots/"
    
    fhr_list = np.arange(first_fhr,last_fhr+fhr_inc,fhr_inc)
    gfs_fhr_list = np.arange(0,246,6)
    gfs_ext_fhr_list = np.arange(252,396,12)
    
    metric_list = ["cenlat","cenlon","ratio","phi","area","t","u"]
    min_data = {}
    max_data = {}
    avg_data = {}
    
    #print "Reading in historical metric data now"
    #for metric in metric_list:
    #    print "Getting "+metric
    #    this_min,this_max,this_avg = get_emetrics_max_min(today,last_fhr,fhr_inc,metric,plot_lev)
    #    min_data[metric] = this_min
    #    max_data[metric] = this_max
    #    avg_data[metric] = this_avg
    print "Opening remote GEFS data "+str(plot_lev)+" and controur level "+str(the_contour)  
    # Just added this if - because its not working for 30 after working for 10 --maybe because we dont need to read in data multiple time
    if plot_lev == 10:

      gefs_files = xr.open_dataset("http://nomads.ncep.noaa.gov:9090/dods/gens/gens"+str(day_str)+"/gep_all_00z")
      g_files = gefs_files["hgtprs"]
      t_files = gefs_files["tmpprs"]
      u_files = gefs_files["ugrdprs"]
      print "Opening remote GFS data at 0.5 deg, 6 hourly res"
      gfs_det_files = xr.open_dataset("http://nomads.ncep.noaa.gov:9090/dods/gfs_0p50/gfs"+str(day_str)+"/gfs_0p50_00z")
      g_files_det = gfs_det_files["hgtprs"]
      t_files_det = gfs_det_files["tmpprs"]
      u_files_det = gfs_det_files["ugrdprs"]
      print "Opening remote GFS data at 1.0 deg, 12 hourly res"
      gfs_ext_files = xr.open_dataset("http://nomads.ncep.noaa.gov:9090/dods/gfs_1p00/gfs"+str(day_str)+"/gfs_1p00_00z")
      g_files_det_ext = gfs_ext_files["hgtprs"]
      t_files_det_ext = gfs_ext_files["tmpprs"]
      u_files_det_ext = gfs_ext_files["ugrdprs"]
    
    ## Select only the data we need for each variable.
    ## GEFS 1.0 deg 0-384 every 6
    print "Selecting and averaging GEFS data."
    print "First, g data"
    g_data = g_files.loc[dict(lat=slice(0,90),lev=plot_lev)]
    print "u data"
    u_data = u_files.loc[dict(lat=65,lev=plot_lev)].mean(dim='lon')
    print "t data"
    t_data = t_files.loc[dict(lat=slice(75,90),lev=plot_lev)]
    t_data = np.average(t_data,axis=2,weights=np.cos(np.radians(t_data['lat'].values)))
    t_data = np.average(t_data,axis=2)
    print "lat and lon"
    lats = g_data["lat"].values
    lons = gefs_files["lon"].values

    ## GFS 1.0 deg data for 252-384 every 12
    print "Selecting and averaging GFS determinsitic data, 12 hourly"
    print "First, g data"
    g_data_det_ext = g_files_det_ext.loc[dict(lat=slice(0,90),lev=plot_lev)]
    print "u data"
    u_data_det_ext = u_files_det_ext.loc[dict(lat=65,lev=plot_lev)].mean(dim='lon')
    print "t data"
    t_data_det_ext = t_files_det_ext.loc[dict(lat=slice(75,90),lev=plot_lev)]
    t_data_det_ext = np.average(t_data_det_ext,axis=1,weights=np.cos(np.radians(t_data_det_ext['lat'].values)))
    t_data_det_ext = np.average(t_data_det_ext,axis=1)
    ## Lat and lon should match, using 1 deg GFS det. data
    
    ## GFS det data for 0-240, 0.5 deg every 6
    print "Selecting and averaging GFS determinsitic data, 6 hourly"
    print "First, g data"
    g_data_det = g_files_det.loc[dict(lat=slice(0,90),lev=plot_lev)]
    print "u data"
    u_data_det = u_files_det.loc[dict(lat=65,lev=plot_lev)].mean(dim='lon')
    print "t data"
    t_data_det = t_files_det.loc[dict(lat=slice(75,90),lev=plot_lev)]
    t_data_det = np.average(t_data_det,axis=1,weights=np.cos(np.radians(t_data_det['lat'].values)))
    t_data_det = np.average(t_data_det,axis=1)
    ## Lat and lon should match, using 1 deg GFS det. data
    lats_0p5 = g_data_det["lat"].values
    lons_0p5 = gfs_det_files["lon"].values
    
    ## Open text file to which to write ellipse metrics.
    fout = open(text_file,'w')
    
    print "Begin plotting forecast maps."
    ## Start making plots
    for fhr in fhr_list:
        emark = []
        eline = []
        cs_temp = []
        mem_list = range(1,21)
        mem_list.append(0)
        mem_list.append(-1)
        print "##########################"
        print "Forecast time: "+str(fhr).zfill(3)
        valid_label = (today+dt.timedelta(hours=fhr)).strftime("Valid: %a %b %d %Y %H%M UTC F"+str(fhr).zfill(3))
        fig = plt.figure(figsize=(12,12),dpi=1200)
        ax = fig.add_subplot(111)
        m = Basemap(projection='ortho',lat_0=90,lon_0=0,resolution='l') 
        m.drawcoastlines(linewidth=0.5, color='#A9A49A')  ##bisque map
        m.drawcountries(linewidth=0.5,color='#A9A49A')
        m.fillcontinents(color='#E5E3E0')  #,lake_color='aqua')  #map fill lighter bisque
        m.drawmeridians(np.arange(0,360,30), color='#B0A4A3')   ##gray latlon
        m.drawparallels(np.arange(-90,90,30),color='#B0A4A3')
        clevs = range(18000,33500,250)
        [x,y] = np.meshgrid(lons,lats)
        [x,y] = m(x,y)
        [x_0p5,y_0p5] = np.meshgrid(lons_0p5,lats_0p5)
        [x_0p5,y_0p5] = m(x_0p5,y_0p5)
        plt.title(str(plot_lev)+"hPa Elliptical Diagnostics\n"+day_label+"\n"+valid_label)
        ## This plots member 0 (control) for now.
        #cs = m.contour(x,y,g_data.isel(ens=0,time=fhr/fhr_inc),levels=clevs,linewidths=1.5,colors='black')
        for mem in mem_list:
            print "********************"
            print "Member: "+str(mem).zfill(2)
            ## If you're looking to plot the deterministic but the GFS doesn't have this FHR, skip this loop
            if mem == -1 and fhr not in np.concatenate((gfs_fhr_list,gfs_ext_fhr_list)):
                print "GFS determinstic data not expected at this forecast hour, skipping..."
                continue
            if mem>0:
                mem_color = "gray"
                mem_lw = 1.5
                mem_ms = 5
                cont_color = "gray"
                cont_lw = 0.5
            elif mem==0:
                print "Control member"
                mem_color = "red"
                mem_lw = 3.0
                mem_ms = 8
                cont_color = "#e59797"
                cont_lw = 1.0
            else: ## mem = -1 (deterministic)
                print "Determinsitic"
                mem_color = "blue"
                mem_lw = 3.0
                mem_ms = 8
                cont_color = "#aec2e2"
                cont_lw = 1.0
            mm = copy.copy(m)
            if mem>=0:
                #g_data_cyc, lons_cyc = addcyclic(g_data.isel(ens=mem,time=fhr/fhr_inc),lons)
                #if mem==1:
                #    [xc,yc]=np.meshgrid(lons_cyc,lats)
                #cs_temp.append(mm.contour(xc,yc,g_data_cyc,levels=clevs,linewidths=cont_lw,colors=cont_color))
                cs_temp.append(mm.contour(x,y,g_data.isel(ens=mem,time=fhr/6),levels=clevs,linewidths=cont_lw,colors=cont_color))
            else:
                if fhr in gfs_fhr_list: ## For Days 0-10, 6 hourly data
                    #g_data_cyc, lons_cyc_0p5 = addcyclic(g_data_det.isel(time=np.where(gfs_fhr_list==fhr)[0][0]),lons_0p5)
                    #[xc_0p5,yc_0p5] = np.meshgrid(lons_cyc_0p5,lats_0p5)
                    #cs_temp.append(mm.contour(xc_0p5,yc_0p5,g_data_cyc,levels=clevs,linewidths=cont_lw,colors=cont_color))
                    cs_temp.append(mm.contour(x_0p5,y_0p5,g_data_det.isel(time=np.where(gfs_fhr_list==fhr)[0][0]),levels=clevs,linewidths=cont_lw,colors=cont_color))
                elif fhr in gfs_ext_fhr_list: ## Days 11-16, 12 hourly data. Else is important here.
                    #g_data_cyc, lons_cyc = addcyclic(g_data_det_ext.isel(time=np.where(gfs_ext_fhr_list==fhr)[0][0]),lons)
                    #cs_temp.append(mm.contour(xc,yc,g_data_cyc,levels=clevs,linewidths=cont_lw,colors=cont_color))
                    cs_temp.append(mm.contour(x,y,g_data_det_ext.isel(time=np.where(gfs_ext_fhr_list==fhr)[0][0]),levels=clevs,linewidths=cont_lw,colors=cont_color))
    
            ## Each contour from clevs has its own array, varying lengths
            ## Find out index of desired contour. Read its segment. Convert back to lat/lon.
            try:
                lev_contour_ind = np.where(np.array(cs_temp[-1].levels)==the_contour)[0][0]
                isoline_list = cs_temp[-1].allsegs[lev_contour_ind]
                print "Number of ellipses: ",len(isoline_list)                                
            except:
                print "Contours not found for",the_contour,"meter level at",plot_lev,"hPa."
                isoline_list = []
            ## Will delete the ensemble's contours, if desired
            if mem > 0: ## If not determinsitic or control:
                for coll in cs_temp[-1].collections:
                    coll.remove()
                
            ## Are there two segments, check if there is a real split
           # if len(isoline_list) > 1:
           #      #find end points of the line
                 

            isocount = 0
            small = 0
            for isoline in isoline_list:
                [iso_lon,iso_lat] = mm(isoline[:,0],isoline[:,1],inverse=True)
                if len(iso_lon)<15:
                    print "-----Not analyzing ellipse with",len(iso_lon),"points, continuing..."  # Check for size!
                    small = 1
                    continue
                #if mem == -1:
                #    print "LONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                #    print iso_lon
                #    print "LATS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                #    print iso_lat

                # Checking to see if contours are closed (0) or need to be joined (1)  Checking before convert lat/lon to radians
                lon_diff = abs(iso_lon[0] - iso_lon[len(iso_lon)-1])
                lat_diff = abs(iso_lat[0] - iso_lat[len(iso_lat)-1])
                join = 0
                if lon_diff > 1 or lat_diff > 1:
                     join = 1
                     print "Diffs lat/lon: ",lat_diff,lon_diff
                iso_lon = np.deg2rad(iso_lon)
                iso_lat = np.deg2rad(iso_lat)
   
        
                ## Now convert it into polar coordinates for ellipse math.
                ex = np.array((np.cos(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))
                ey = np.array((np.sin(iso_lon)*np.cos(iso_lat))/(1+np.sin(iso_lat)))

                ## does this contour include the pole?
                overpole = point_inside_polygon(0,0,ex,ey)  #returns true if poly includes the pole, false if not 
                print "It is",overpole, "that the contour includes the pole"
                   
                ### contours break if they don't include the pole and need to be reconnected here...
                ### first check that contour doesn't include pole, there is more than one contour,... 
                ###        ...the contours are not a real split vortex and we didn't get rid of one because it was small...
                if not overpole and len(isoline_list) > 1 and join > 0 and small < 1:
                     ex2 = ex
                     ey2 = ey
                     if isocount > 0:
                          ex = np.append(ex,ex2)
                          ey = np.append(ey,ey2)
                          print "Point didn't include pole - add to these to previous set"
                     else:
                          print "Point didn't include pole - keeping these to add to next set"
                          isocount = 1
                          continue

                print "Running ellipse diagnostic now"
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
                print "Center of ellipse:",cenlat,"N",cenlon,"E"
        
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
        
                ## Plot full ellipse and center mark #############
                exm,eym = m(elons,elats)
                exy = zip(exm,eym)
                eline.append(Polygon(exy,closed=True,color=mem_color,fill=False,lw=mem_lw))
                ax.add_patch(eline[-1])
                mcx,mcy = m(cenlon,cenlat)
                emark.append(m.plot(mcx,mcy,marker='o',markerfacecolor=mem_color,markersize=mem_ms))
        
                ephi = np.rad2deg(ephi)
                ## If a1 shorter than b1, adjust phi to correct the orientation angle
                if a1gc < b1gc:
                    ephi -= 90
                print "Emetrics phi:",ephi
                if ephi < -45:
                    ephi += 180
                
                ratio = a1gc/b1gc
                if ratio < 1.0:
                    ratio = 1.0/ratio
        
                size = math.pi*a1gc*b1gc
        
                ## Print the text file, format:
                ##  hhstr+" "+ddstr+" "+mmstr+" "+yystr+" "+fff+" "+mbr+" "+cenlondeg+" "+cenlatdeg+" "+a1gcircle+" "+b1gcircle+" "+phideg+" "+ufcast+" "+tfcast
                ## 00 21 Oct 2016 000 0 -159.609 88.6679 1738.34 2479.64 -89.5687 26.7172 205.673
                fhr_idx = np.where(fhr_list==fhr)[0]
                mem_idx = mem #np.where(u_data["ens"].values==mem)[0]
                if mem > -1:
                    fhr_time = today+dt.timedelta(hours=fhr)
                    fhr64 = np.datetime64(fhr_time)
                    u_data_mem = u_data.loc[dict(time=fhr64,ens=mem+1)].values
                    t_data_mem = t_data[mem_idx,fhr_idx][0]
                    mem_str = str(mem).zfill(2)
                else:
                    mem_str = 99
                    if fhr<=240:
                        u_data_mem = u_data_det[fhr_idx].values[0]
                        t_data_mem = t_data_det[fhr_idx][0]
                    else:
                        ext_fhr_idx = np.where(gfs_ext_fhr_list==fhr)[0]
                        u_data_mem = u_data_det_ext[ext_fhr_idx].values[0]
                        t_data_mem = t_data_det_ext[ext_fhr_idx][0]
                estr = today.strftime("%H %Y%m%d")+" "+str(fhr).zfill(3)+" "+str(mem_str)+" "+str(cenlon)+" "+str(cenlat)+" "+str(a1gc)+" "+str(b1gc)+" "+str(ephi)+" "+str(u_data_mem)+" "+str(t_data_mem)
                fout.write(estr+"\n")
    
        ## Setting up contour legend.
        labels = ['GFS Heights','GEFS Ctrl','GEFS Mem','GFS Ellipse','GEFS Ctrl','GEFS Mem','GFS center','GEFS Ctrl center','GEFS Mem center']
        if fhr in np.concatenate((gfs_fhr_list,gfs_ext_fhr_list)):
            cs_temp[-1].collections[0].set_label(labels[0])
            eline[-1].set_label(labels[3])
            cs_temp[-2].collections[0].set_label(labels[1])
            eline[-2].set_label(labels[4])
        else:
            cs_temp[-1].collections[0].set_label(labels[1])
            eline[-1].set_label(labels[4])
        cs_temp[0].collections[0].set_label(labels[2])
        eline[0].set_label(labels[5])
        #emark[-1].collections[0].set_label(labels[6])
        #emark[-2].collections[0].set_label(labels[7])
        #emark[0].collections[0].set_label(labels[8])
        handles,leg_labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1],leg_labels[::-1],loc='lower right',ncol=1,markerscale=1.8)
    
        ## Save figures
    #    plt.savefig("test_ellipse_all_F"+str(fhr).zfill(3)+".png",format='png')
        plt.savefig(plots_loc+"ellipse"+str(plot_lev)+"_"+str(fhr/6)+".png",format='png')
        plt.close()
    

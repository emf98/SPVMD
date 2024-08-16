#!/usr/bin/python
# -*- coding: utf-8 -*-
#Solve for the best fit for an ellipse then plot it!
# Generated: 7 Aug 2015 by A Lang
# Most Recent: 25 July 2024 E Fernandez

# Uses x.txt and y.txt produced by NCL to find lat/lon of 30000m 10-hPa or 23000m 30-hPa contour which 
# then converts lat/lon to cartesian coords with N.Pole at origin
# This python script finds ellipse in cart. coords, then writes cart. coord of ellipse as xx.txt and yy.txt to be used by NCL
# 

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

import os
import numpy as np
from numpy.linalg import eig, inv


def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1

#The code contains the original and the updated ways of solving for the ellipse
#Old way ---->    Authors note it is numerically unstable if not implemented in matlab
#    E, V =  eig(np.dot(inv(S), C))
#    print("E: ", E
#    print("V: ",V  
#    n = np.argmax(np.abs(E))
#    print("n: ",n
#    a = V[:,n]
   
#New way ----> numerically stabler in C [gevec, geval] = eig(S,C);
#    print("new method"
# First break matrix into blocks
   
    tmpA = S[0:3,0:3]
    tmpB = S[0:3,3:6]
    tmpC = S[3:,3:]
   
    tmpD = C[0:3,0:3]
    tmpE = np.dot(inv(tmpC), tmpB.conj().T)
    tmpF = np.dot(tmpB,tmpE)
    eval_x, evec_x = eig(np.dot(inv(tmpD), (tmpA - tmpF)))
  
# Find the negative (as det(tmpD) < 0) eigenvalue
    I = np.argmax(eval_x)
  
#Extract eigenvector corresponding to negative eigenvalue
    a = evec_x[:,I]
  
# Recover the bottom half...
    evec_y = np.dot(-1*tmpE, a)
    a = np.concatenate((a,evec_y))  

      
    return a


# These three functions are from the old way, using the quadratic method 
# of solving for the ellipse 
# -> ellipse_axis_length errors out when down1/2 contain the sqrt of a (-)
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

# This function is the updated way of solving for an ellipse (geometric not quadratic)
# Works in cases where down1/2 in above way is negative
def get_ellipse_metrics( a ):
    thtarad = .5 * np.arctan2(a[1], a[0]-a[2])
    cost = np.cos(thtarad)
    sint = np.sin(thtarad)
    sinsq = sint*sint
    cossq = cost*cost
    cossin = sint*cost
    
    Ao = a[5]
    Au = a[3]*cost + a[4]*sint
    Av = -1*a[3]*sint + a[4]*cost
    Auu = a[0]*cossq + a[2]*sinsq + a[1]*cossin
    Avv = a[0]*sinsq + a[2]*cossq - a[1]*cossin
    
    tuCenter = -1*Au/(2*Auu)
    tvCenter = -1*Av/(2*Avv)
    wCenter = Ao - Auu*tuCenter*tuCenter - Avv*tvCenter*tvCenter
    uCenter = tuCenter*cost - tvCenter*sint
    vCenter = tuCenter*sint + tvCenter*cost
    
    Ru = -1*wCenter/Auu
    Rv = -1*wCenter/Avv
    
    Ru = np.sqrt(np.abs(Ru))*np.sign(Ru)  
    Rv = np.sqrt(np.abs(Rv))*np.sign(Rv)
    
    return np.array([uCenter, vCenter, Ru, Rv, thtarad])



def get_numbers_from_file(file_name):
    file = open(file_name, "r")
#    strnumbers = file.read().split()
#    return map(float, strnumbers)
    return np.loadtxt(file)

def get_ellipse_coords(a, b, x, y, angle, k=2):
    # Source: scipy-central.org/item/23/1/plot-an-ellipse
    # License: Creative Commons Zero (almost public domain) http://scpyce.org/cc0
    # Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    # given at http://en.wikipedia.org/wiki/Ellipse
    # k = 1 (2) means 361 (721) points (degree by degree)
    # a = major axis distance,
    # b = minor axis distance,
    # x = offset along the x-axis
    # y = offset along the y-axis
    # angle = clockwise rotation [in degrees] of the ellipse;
    #    * angle=0  : the ellipse is aligned with the positive x-axis
    #    * angle=30 : rotated 30 degrees clockwise from positive x-axis
    
    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])
 
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    
    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts


####------------------------------------------------------------------------
####---------- Start Code --------------------------------------------------
#### Define functions first:
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

def fitEllipseContour(x,y):

    a = fitEllipse(x,y)
    #uses old way
    center = ellipse_center(a)
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    #print("old center = "+str(center))
    #print("old angle of rotation = "+str(phi*180/np.pi))
    #print("old axes = "+ str(axes))
    
    #uses new way to get ellipse metrics
    edata = get_ellipse_metrics( a )
    center = np.array([edata[0],edata[1]])
    phi = edata[4]
    axes = np.array([edata[2],edata[3]])
    
    phideg = -phi * 180.0/np.pi
    #print("center = "+str( center))
    #print("angle of rotation = "+str(phi*180/np.pi))
    #print("axes = "+str(axes))
    
    a, b = axes

    axes_center_phi = str(a)+" "+str(b)+" "+str(center[0])+" "+str(center[1])+" "+str(phi)
    #print(axes_center_phi)
    pts = get_ellipse_coords(a,b,center[0],center[1],phideg,2)
    circ = get_ellipse_coords(.5,.5,0,0,0,2)   # plot a circle equivalent to the 45N circle
    
    xx = pts[:,0]
    yy = pts[:,1]
    
    return xx,yy,a,b,center[0],center[1],phi

    ##xx_ret = map(str,xy) for xy in (xx,xx)
    ##yy_ret =

    #np.savetext('/al11/andrea/research/Elliptical/xxyy.txt', (xx,yy), fmt="%d")
    #with open (dirpath + '/xx.txt','w') as fx:
    #    fx.writelines("\n".join(map(str,xy)) for xy in (xx,xx)) 
    #
    #with open (dirpath + '/yy.txt','w') as fy:
    #    fy.writelines("\n".join(map(str,yx)) for yx in (yy,yy)) 
    
    #from pylab import *
    #plot(x,y, color = 'blue', linewidth=2)    #plot the geopotential contour that was provided converted to cartesian coords
    #plot(center[0],center[1], 'ro', markersize = 8)  # plot the center of that closed contour and ellipse as a 'o'
    #plot(xx,yy, color = 'red', linewidth=3) # plot the best fit ellipse for that contour in red
    #plot(circ[:,0],circ[:,1],'green')  # plot a circle equivalent to the 45N circle
    #plot(0,0,'g+',markersize=25) #plot the north poles
    #show()



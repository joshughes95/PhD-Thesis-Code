# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:51:52 2024

@author: tscw79
"""

from os.path import dirname, join as pjoin
import scipy.io as sio
import scipy
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askdirectory
import csv
import os
import pandas as pd
from skimage import io
import cv2
from itertools import zip_longest
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm

#########################################################

Input_directory = askdirectory(title='Select input Directory')

Results_directory = askdirectory(title='Select Results Directory')

filename = 'LightOx78 500uM CCl4 pump 100 nJ at 400nm probe 2200cm-1 1574cm-1 Flow'

#########################################################

def Lifetime(x, A, B, C):
    y = A*np.exp(-B*x)+C
    return y

def Gauss(x, A, B, C):
    y = A*np.exp(-((x-B)**2)/(C**2))
    return y

def GaussPlus(x, A, B, C, D):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D
    return y
    
def DoubleGauss(x, A, B, C, D, E, F):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2))
    return y

def DoubleGaussPlus(x, A, B, C, D, E, F, G):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2)) + G
    return y

def TripleGauss(x, A, B, C, D, E, F, G, H, I):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2)) + G*np.exp(-((x-H)**2)/(I**2))
    return y

def TripleGaussPlus(x, A, B, C, D, E, F, G, H, I, K):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2)) + G*np.exp(-((x-H)**2)/(I**2)) + K
    return y

def QuadGauss(x, A, B, C, D, E, F, G, H, I, J, K, L):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2)) + G*np.exp(-((x-H)**2)/(I**2)) + J*np.exp(-((x-K)**2)/(L**2))
    return y

def QuadGaussPlus(x, A, B, C, D, E, F, G, H, I, J, K, L, M):
    y = A*np.exp(-((x-B)**2)/(C**2)) + D*np.exp(-((x-E)**2)/(F**2)) + G*np.exp(-((x-H)**2)/(I**2)) + J*np.exp(-((x-K)**2)/(L**2)) + M
    return y

#########################################################
df_list = []
sheet_names = []
filenames = []

df = pd.read_csv(Input_directory+'\\'+filename+'.csv')

filenames.append(filename)
wavenumbers = df.iloc[:,1].to_numpy().astype(float).transpose()
pixels = df.iloc[:,0].to_numpy().astype(int).transpose()

time = df.columns[2:66]
heat_map = df.iloc[:,2:66].transpose()
heat_map.replace('#NAME?', float(np.nan), inplace=True)
heat_map.replace(np.inf, float(np.nan), inplace=True)
heat_map = heat_map.astype(float)

aspect_ratio = len(wavenumbers)/len(time)

########################### show Heat Map ################      
plt.imshow(heat_map, cmap='hot', interpolation='none', vmin= -0.00075, vmax=0.00075)
plt.title(filename)
#plt.xticks(wavenumbers, wavenumbers)
#plt.yticks(np.arange(0,len(time),1), time)
plt.gca().set_aspect(aspect_ratio)
plt.colorbar()
#plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'_heat_map')
plt.show()

########### FT-IR to show peaks vs wavenumbers ###########
# x = np.r_[wavenumbers[1]:wavenumbers[128]:(wavenumbers[129]-wavenumbers[1])/128
#           , wavenumbers[129]:wavenumbers[256]:(wavenumbers[256]-wavenumbers[129])/128]

fig,(ax,ax2) = plt.subplots(1, 2, sharey=True)            
# plot the same data on both axes
scale_factor = max(heat_map.iloc[20])/max(heat_map.iloc[55])
ax.plot(wavenumbers, heat_map.iloc[20])
ax2.plot(wavenumbers, heat_map.iloc[20], label='time='+str(time[20])+'ps')
ax.plot(wavenumbers, heat_map.iloc[55])
ax2.plot(wavenumbers, heat_map.iloc[55], label='time='+str(time[55])+'ps')

# zoom-in / limit the view to different portions of the data
ax.set_xlim(wavenumbers[0], wavenumbers[127])
ax2.set_xlim(wavenumbers[128], wavenumbers[255])

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
#ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_right()
plt.xlabel('Wavenumbers cm-1')
ax.set_ylabel(r'$\Delta$ Absorbance')
ax.set_ylim(-0.0002, 0.0004)
ax2.set_ylim(-0.0002, 0.0004)
plt.title(filename)
plt.legend()
#plt.savefig(Results_directory+'\\'+compound+'\\'+solvent + '\\'+str(filename)+'_FT_IR')
plt.show()

######################## Select bleach points to fit a Gaussian ###############

pixel_bleach = [34, 39, 77, 84, 185]
pixel_s1 = [21, 79]
pixel_t1 = [35, 66, 100]

ymin = -0.0002
ymax = 0.0004

fig,(ax,ax2) = plt.subplots(1, 2, sharey=True)            
# plot the same data on both axes
scale_factor = max(heat_map.iloc[20])/max(heat_map.iloc[55])
ax.plot(wavenumbers, heat_map.iloc[20])
ax2.plot(wavenumbers, heat_map.iloc[20], label='time='+str(time[20])+'ps')
ax.plot(wavenumbers, heat_map.iloc[55])
ax2.plot(wavenumbers, heat_map.iloc[55], label='time='+str(time[55])+'ps')

# zoom-in / limit the view to different portions of the data
ax.set_xlim(wavenumbers[0], wavenumbers[127])
ax2.set_xlim(wavenumbers[128], wavenumbers[255])

ax.vlines(x=wavenumbers[pixel_bleach], ymin=ymin, ymax=ymax, color ='purple', linestyle='--')
ax2.vlines(x=wavenumbers[pixel_bleach], ymin=ymin, ymax=ymax, color= 'purple', linestyle='--')

ax.vlines(x=wavenumbers[pixel_s1], ymin=ymin, ymax=ymax, color ='green', linestyle='--')
ax2.vlines(x=wavenumbers[pixel_s1], ymin=ymin, ymax=ymax, color= 'green', linestyle='--')

ax.vlines(x=wavenumbers[pixel_t1], ymin=ymin, ymax=ymax, color ='pink', linestyle='--')
ax2.vlines(x=wavenumbers[pixel_t1], ymin=ymin, ymax=ymax, color= 'pink', linestyle='--')

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
#ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_right()
plt.xlabel('Wavenumbers cm-1')
ax.set_ylabel(r'$\Delta$ Absorbance')
ax.set_ylim(ymin, ymax)
ax2.set_ylim(ymin, ymax)
plt.title(filename)
plt.legend()
#plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'_FT_IR')
plt.show()

########################### Fit Gaussians ####################################

# Bleach in Triplet states
s1_average = np.nanmean(heat_map.iloc[16:20], axis=0)
t1_average = np.nanmean(heat_map.iloc[54:64], axis=0)

bleach1 = pixel_bleach[2] # select bleach
bleach2 = pixel_bleach[3] # select bleach
s1 = pixel_s1[1]
s2 = pixel_s1[1]
t1 = pixel_t1[1]
t2 = pixel_t1[2]

xdata = pixels[15:55]
ydata = t1_average[15:55]

plt.plot(xdata, ydata) # visualise it
plt.show()

#####IF CLOSE TO ANOTHER BAND, USE DOUBLE GAUSS
# SINGLE GAUSS
bleach1 = 85
boundary = ([-0.001, bleach1-5, 0], [0.000, bleach1+5, 50]) # boundaries and guesses
p0 = np.array([-0.0001, bleach1, 5])


parameters, covariance = curve_fit(Gauss, xdata, ydata, bounds=boundary, p0=p0)
fit = Gauss(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[53]+'-'+time[63]+'ps')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, 
                                      parameters[0], parameters[1], parameters[2]), linestyle='--', label='fit')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - Gauss(xdata,
                            parameters[0], parameters[1], parameters[2]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

# SINGLE GAUSS PLUS
xdata = pixels[65:100]
ydata = t1_average[65:100]

bleach1=232.75
s1=200

boundary = ([-0.001, bleach1-0.16, 0, -0.0001], [0.000, bleach1+0.16, 50, 0.0001]) # boundaries and guesses
p0 = np.array([-0.00003, bleach1, 3, 0])


parameters, covariance = curve_fit(GaussPlus, xdata, ydata, bounds=boundary, p0=p0)
fit = GaussPlus(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[53]+'-'+time[63]+'ps')
axs[0].plot(wavenumbers[xdata], GaussPlus(xdata, 
                                      parameters[0], parameters[1], parameters[2], parameters[3]), label='fit')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - GaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2], parameters[3]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

#DOUBLE GAUSS
xdata = pixels[15:55]
ydata = s1_average[15:55]

bleach1=80
s1=70

boundary = ([-0.001, s1-5, 0, -0.001, bleach1-5, 0], 
            [0.001, s1+5, 50, 0.000, bleach1+5, 50]) # boundaries and guesses
p0 = np.array([0.0004, s1, 5, -0.0002, bleach1, 5])


parameters, covariance = curve_fit(DoubleGauss, xdata, ydata, bounds=boundary, p0=p0)
fit = DoubleGauss(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[16]+'-'+time[20]+'ps')
axs[0].plot(wavenumbers[xdata], DoubleGauss(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5]), label='fit')

axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[0], parameters[1], parameters[2]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]), linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - DoubleGauss(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

#DOUBLE GAUSS + C
xdata = pixels[60:90]
ydata = t1_average[60:90]
bleach1 = 79.8
boundary = ([-0.001, t1-5, 0, -0.001, bleach1-0.26, 0, -0.0001], 
            [0.001, t1+5, 20, 0.000, bleach1+0.26, 50, 0.0001]) # boundaries and guesses
p0 = np.array([0.0003, t1, 10, -0.0001, bleach1, 4, 0])


parameters, covariance = curve_fit(DoubleGaussPlus, xdata, ydata, bounds=boundary, p0=p0)
fit = DoubleGaussPlus(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[50]+'-'+time[60]+'ps')
axs[0].plot(wavenumbers[xdata], DoubleGaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5], parameters[6]), label='fit')

axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[0], parameters[1], parameters[2]) + parameters[6], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]) + parameters[6], linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - DoubleGaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5], parameters[6]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

#TRIPLE GAUSS
xdata = pixels[60:110]
ydata = s1_average[60:110]

s1=82
boundary = ([-0.001, bleach1-5, 0, -0.001, s1-5, 0, -0.001, bleach2-5, 0],
            [0.001, bleach1+5, 50, 0.001, s1+5, 100, 0.001, bleach2+5, 50]) # boundaries and guesses
p0 = np.array([-0.00025, bleach1, 5, 0.0003, s1, 20, -0.0002, bleach2, 5])


parameters, covariance = curve_fit(TripleGauss, xdata, ydata, bounds=boundary, p0=p0)
fit = TripleGauss(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[16]+'-'+time[20]+'ps')
axs[0].plot(wavenumbers[xdata], TripleGauss(xdata,parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8]))

axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[0], parameters[1], parameters[2]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[6], parameters[7], parameters[8]), linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - TripleGauss(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

#TRIPLE GAUSS PLUS
xdata = pixels[60:105]
ydata = s1_average[60:105]

bleach1 = 30.38
bleach2 = 38.66
boundary = ([-0.001, s1-5, 0, -0.001, bleach1-5, 0, -0.001, bleach2-5, 0, -0.0001],
            [0.001, s1+5, 50, 0.001, bleach1+5, 100, 0.001, bleach2+5, 50, 0.0001]) # boundaries and guesses
p0 = np.array([0.0003, s1, 5, -0.0002, bleach1, 2, -0.0002, bleach2, 3, 0.0001])


parameters, covariance = curve_fit(TripleGaussPlus, xdata, ydata, bounds=boundary, p0=p0)
fit = TripleGaussPlus(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[16]+'-'+time[20]+'ps')
axs[0].plot(wavenumbers[xdata], TripleGaussPlus(xdata,parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8], parameters[9]))

axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[0], parameters[1], parameters[2]) + parameters[9], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]) + parameters[9], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[6], parameters[7], parameters[8]) + parameters[9], linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - TripleGaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8], parameters[9]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()


#QUAD GAUSS
xdata = pixels[66:110]
ydata = t1_average[66:110]

boundary = ([-0.001, t1-5, 0, -0.0005, bleach1-5, 0, -0.001, bleach2-5, 0, -0.001, t2-5, 0],
            [0.001, t1+5, 50, -0.00001, bleach1+5, 50, 0.001, bleach2+5, 50, 0.001, t2+5, 50]) # boundaries and guesses
p0 = np.array([0.0001, t1, 5, -0.00005, bleach1, 3, -0.0001, bleach2, 5, 0.00005, t2, 5])


parameters, covariance = curve_fit(QuadGauss, xdata, ydata, bounds=boundary, p0=p0)
fit = QuadGauss(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[50]+'-'+time[60]+'ps')
axs[0].plot(wavenumbers[xdata], QuadGauss(xdata,parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11]))

axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[0], parameters[1], parameters[2]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[6], parameters[7], parameters[8]), linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[9], parameters[10], parameters[11]), linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - QuadGauss(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

#QUAD GAUSS PLUS
xdata = pixels[68:105]
ydata = t1_average[68:105]

bleach1 = 34.69
bleach2 = 41.96
boundary = ([-0.001, t1-5, 0, -0.0005, bleach1-5, 0, -0.001, bleach2-5, 0, -0.001, t2-5, 0, -0.001],
            [0.001, t1+5, 50, -0.00001, bleach1+5, 50, 0.001, bleach2+5, 50, 0.001, t2+5, 50, 0.001]) # boundaries and guesses
p0 = np.array([0.0003, t1, 1, -0.0005, bleach1, 1, -0.0001, bleach2, 5, 0.0005, t2, 5, 0.0001])


parameters, covariance = curve_fit(QuadGaussPlus, xdata, ydata, bounds=boundary, p0=p0)
fit = QuadGaussPlus(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[53]+'-'+time[63]+'ps')
axs[0].plot(wavenumbers[xdata], QuadGaussPlus(xdata,parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11], parameters[12]))

axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[0], parameters[1], parameters[2]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[6], parameters[7], parameters[8]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[9], parameters[10], parameters[11]) + parameters[12], linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - QuadGaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11], parameters[12]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()

xdata = pixels[60:110]
ydata = s1_average[60:110]

s1 = 77
s2 = 85
boundary = ([-0.001, s1-5, 0, -0.0005, bleach1-5, 0, -0.001, bleach2-5, 0, -0.001, s2-5, 0, -0.001],
            [0.001, s1+5, 50, -0.00001, bleach1+5, 50, 0.001, bleach2+5, 50, 0.001, s2+5, 50, 0.001]) # boundaries and guesses
p0 = np.array([0.0003, s1, 5, -0.0002, bleach1, 5, -0.0002, bleach2, 5, 0.0003, s2, 5, 0.0001])


parameters, covariance = curve_fit(QuadGaussPlus, xdata, ydata, bounds=boundary, p0=p0)
fit = QuadGaussPlus(xdata, *parameters)
r_sq = r2_score(ydata, fit)
errors = np.sqrt(np.diag(covariance))

fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
axs[0].plot(wavenumbers[xdata], ydata, label='Time Delay = '+time[16]+'-'+time[20]+'ps')
axs[0].plot(wavenumbers[xdata], QuadGaussPlus(xdata,parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11], parameters[12]))

axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[0], parameters[1], parameters[2]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata, parameters[3], parameters[4], parameters[5]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[6], parameters[7], parameters[8]) + parameters[12], linestyle='--')
axs[0].plot(wavenumbers[xdata], Gauss(xdata,parameters[9], parameters[10], parameters[11]) + parameters[12], linestyle='--')
axs[0].set_ylabel(r'$\Delta$ Absorbance')

axs[1].plot(wavenumbers[xdata], ydata - QuadGaussPlus(xdata,
                            parameters[0], parameters[1], parameters[2],
                            parameters[3], parameters[4], parameters[5],
                            parameters[6], parameters[7], parameters[8],
                            parameters[9], parameters[10], parameters[11], parameters[12]))
axs[1].set_ylabel('Residual')
axs[1].set_xlabel('Wavenumbers cm$^{-1}$')
fig.legend()
plt.show()






######################## Find maxima and plot peak wavenumber vs time #########
#peaks = scipy.signal.find_peaks(heat_map.iloc[20][128:256])
max_peak = np.where(max(heat_map.iloc[20][128:256]) == heat_map.iloc[20])[0][0]
max_peak_time = np.where(max(heat_map.iloc[:][max_peak]) == heat_map.iloc[:][max_peak])[0][0]

plt.plot(time.to_numpy().astype(float), heat_map.iloc[:][max_peak])
plt.xlabel('Time Delay / ps')
plt.ylabel(r'$\Delta$ Absorbance')
plt.title(filename)
#fig.set_major_locator(plt.MaxNLocator(10))
plt.locator_params(axis='x',nbins=10)
plt.show()

# Curve fit to extract lifetime
xdata = time.to_numpy()[max_peak_time:].astype(float)
ydata = heat_map.iloc[:][max_peak].to_numpy().transpose()[max_peak_time:]
Lifetime_p0 = np.array([heat_map.iloc[20][max_peak], 0.0005, 0])
parameters, covariance = curve_fit(Lifetime, xdata, ydata, p0 = Lifetime_p0, maxfev=2000)
y_pred = Lifetime(xdata, *parameters)
r_sq = r2_score(ydata, y_pred)
errors = np.sqrt(np.diag(covariance))

plt.plot(time.to_numpy().astype(float), heat_map.iloc[:][max_peak], label='Measured @ '+str(round(wavenumbers[max_peak]))+'cm-1')
plt.plot(xdata, y_pred, label=r'Best Fit, $\tau$='+ str(round(1/parameters[1]))+'ps')
plt.xlabel('Time Delay / ps')
plt.ylabel(r'$\Delta$ Absorbance')
plt.title(filename)
#fig.set_major_locator(plt.MaxNLocator(10))
plt.locator_params(axis='x',nbins=10)
plt.legend()
plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'liftime')
plt.show()

plt.plot(time, heat_map.iloc[:][max_peak], label='Measured @ '+str(round(wavenumbers[max_peak]))+'cm-1')
y_pred = Lifetime(time.to_numpy().astype(float), *parameters)
plt.plot(time, y_pred, label=r'Best Fit, $\tau$='+ str(round(1/parameters[1]))+'ps')
plt.xlabel('Time Delay / ps')
plt.ylabel(r'$\Delta$ Absorbance')
plt.title(filename)
#fig.set_major_locator(plt.MaxNLocator(10))
plt.locator_params(axis='x',nbins=10)
plt.legend()
plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'lifetime_difxscale')
plt.show()

fit_A.append(parameters[0])
fit_B.append(parameters[1])
fit_C.append(parameters[2])

fit_A_err.append(errors[0])
fit_B_err.append(errors[1])
fit_C_err.append(errors[2])

r_sqs.append(r_sq)
lifetime.append(1/parameters[1])
lifetime_err.append(np.sqrt((1/parameters[1])**2 - (1/(parameters[1] + errors[1]))**2))

    df_lifetime = pd.DataFrame({'filename': filenames, 'fit_A':fit_A, 'fit_A_err':fit_A_err,
                                'fit_B':fit_B, 'fit_B_err':fit_B_err,
                                'fit_C':fit_C, 'fit_C_err':fit_C_err,'r_sq': r_sqs,
                                'lifetime':lifetime, 'lifetime_err':lifetime_err})
    
    df_list.append(df_lifetime)
    
Excelwriter = pd.ExcelWriter(Results_directory+'\\Singlet Lifetimes.xlsx',engine="xlsxwriter")
for i, df in enumerate (df_list):    
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
    
Excelwriter.save()

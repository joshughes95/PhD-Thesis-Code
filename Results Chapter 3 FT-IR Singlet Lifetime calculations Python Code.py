# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:57:39 2023

@author: tscw79
"""


from os.path import dirname, join as pjoin
import scipy.io as sio
import scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import sys
from tkinter import Tk
from tkinter.filedialog import askdirectory
import csv
import os
import pandas as pd
from skimage import io
#from moviepy.editor import VideoClip
#from moviepy.video.io.bindings import mplfig_to_npimage
import cv2
from cv2_rolling_ball import subtract_background_rolling_ball
from itertools import zip_longest
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
from lmfit import Parameters, minimize, report_fit
from matplotlib.widgets import SpanSelector
import symfit as sf
from symfit import parameters, variables, Fit, exp

###############################################################################

Input_directory = askdirectory(title='Select input Directory')

Results_directory = askdirectory(title='Select Results Directory')

###############################################################################

def Lifetime(x, A, B, C):
    y = A*np.exp(-B*x)+C
    return y

def LifetimeMonoExp_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    A = params[f'A_{i+1}']
    #B = params[f'B_{i+1}']
    B = params['B']
    C = params[f'C_{i+1}']

    return Lifetime(x, A, B, C)


def objectiveMono(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - LifetimeMonoExp_dataset(params, i, x[i])

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

def LifetimeDoubleExp(x, A, B, C, D, E):
    y = A*np.exp(-B*x)+C*np.exp(-D*x)+E
    return y

def LifetimeDoubleExp_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    A = params[f'A_{i+1}']
    B = params['B']
    C = params[f'C_{i+1}']
    D = params['D']
    E = params[f'E_{i+1}']
    return LifetimeDoubleExp(x, A, B, C, D, E)


def objectiveDouble(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - LifetimeDoubleExp_dataset(params, i, x[i])

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

def LifetimeTripleExp(x, A, B, C, D, E, F, G):
    y = A*np.exp(-B*x)+C*np.exp(-D*x)+E*np.exp(-F*x)+G
    return y

def LifetimeTripleExp_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    A = params[f'A_{i+1}']
    B = params['B']
    C = params[f'C_{i+1}']
    D = params['D']
    E = params[f'E_{i+1}']
    F = params['F']
    G = params[f'G_{i+1}']
    return LifetimeTripleExp(x, A, B, C, D, E, F, G)


def objectiveTriple(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - LifetimeTripleExp_dataset(params, i, x[i])

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

###############################################################################
pixels = np.arange(0,256,1)

for compound in os.listdir(Input_directory):
    df_list = []
    sheet_names = []
    for solvent in os.listdir(Input_directory+'\\'+compound):
        sheet_names.append(solvent)
        

        
        filenames = []
        filenames_4_df = []
        peak_4_df = []
        peak_4_df_pixel = []
        
        fit_A_trip = []
        fit_B_trip = []
        fit_C_trip = []
        fit_D_trip = []
        fit_E_trip = []
        fit_F_trip = []
        fit_G_trip = []
        fit_A_trip_err = []
        fit_B_trip_err = []
        fit_C_trip_err = []
        fit_D_trip_err = []
        fit_E_trip_err = []
        fit_F_trip_err = []
        fit_G_trip_err = []
        residuals_trip = []
        chi_sq_trip =[]
        lifetime1_trip = []
        lifetime1_trip_err = []
        lifetime2_trip = []
        lifetime2_trip_err = []
        lifetime3_trip = []
        lifetime3_trip_err = []
        av_amp_lifetime_trip = []
        av_int_lifetime_trip = []
        
        fit_A = []
        fit_B = []
        fit_C = []
        fit_D = []
        fit_E = []
        fit_A_err = []
        fit_B_err = []
        fit_C_err = []
        fit_D_err = []
        fit_E_err = []
        residuals = []
        chi_sq =[]
        lifetime1 = []
        lifetime1_err = []
        lifetime2 = []
        lifetime2_err = []
        av_amp_lifetime = []
        av_int_lifetime = []
        
        fit_A_mono = []
        fit_B_mono = []
        fit_C_mono = []
        fit_A_mono_err = []
        fit_B_mono_err = []
        fit_C_mono_err = []
        residuals_mono = []
        chi_sq_mono = []
        lifetime_mono = []     
        lifetime_mono_err = []
        
        for filename in os.listdir(Input_directory+'\\'+compound+'\\'+solvent):
            if 'EO' not in filename:
                df_output = pd.DataFrame(columns=('Concentration', 'Abs')) # df for each conc
                        
                df = pd.read_csv(Input_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename))
                filename = os.path.splitext(filename)[0]
                filenames.append(filename)
                wavenumbers = df.iloc[:,1].to_numpy().astype(float).transpose()
                time = df.columns[2:].to_numpy()
                time[np.where(time == '0.6')[0]] = 0.0
                heat_map = df.iloc[:,2:].transpose()
                heat_map.replace('#NAME?', float(np.nan), inplace=True)
                heat_map.replace(np.inf, float(np.nan), inplace=True)
                heat_map = heat_map.astype(float)
                
                aspect_ratio = len(wavenumbers)/64
        
    ########################### show Heat Map #####################################        
                plt.imshow(heat_map[:64], cmap='hot', interpolation='none', aspect='auto')
                plt.title(filename)
                plt.xticks(np.arange(28,268,40), [round(y) for y in (wavenumbers[np.arange(28,268,40)])])
                plt.yticks(np.arange(2,72,10), [round(y) for y in ([round(y) for y in ((time[np.arange(2,72,10)]).astype(float))])])
                plt.gca().set_aspect(aspect_ratio)
                plt.colorbar(label=r'$\Delta$ Absorbance')
                plt.xlabel('Wavenumbers / cm$^{-1}$')
                plt.ylabel('Time Delay / ps')
                plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'_heat_map', bbox_inches='tight')
                plt.show()
                
    ########################## IR to show peaks vs wavenumbers #################
    
                fig,(ax,ax2) = plt.subplots(1, 2, sharey=False)            
                # plot the same data on both axes
                
                time_delays = np.array([12,15,20,25,30,35,55])
                
                for td in time_delays:
                    ax.plot(wavenumbers[20:107], heat_map.iloc[td][20:107])
                    ax2.plot(wavenumbers[128:255], heat_map.iloc[td][128:255], label='time='+str(time[td])+'ps')
    
                
                # zoom-in / limit the view to different portions of the data
                #ax.set_xlim(wavenumbers[0], wavenumbers[127])
                #ax2.set_xlim(wavenumbers[128], wavenumbers[255])
                
                # hide the spines between ax and ax2
                #ax.spines['right'].set_visible(False)
                #ax2.spines['left'].set_visible(False)
                ax.yaxis.tick_left()
                #ax.tick_params(labeltop='off') # don't put tick labels at the top
                ax2.yaxis.tick_right()
                ax.set_xlabel('Wavenumbers / cm$^{-1}$')
                ax2.set_xlabel('Wavenumbers / cm$^{-1}$')
                ax.set_ylabel(r'$\Delta$ Absorbance')
                plt.title(filename)
                plt.legend(loc='center right', bbox_to_anchor=(2.2, 0.5),
                           ncol=1, fancybox=True, shadow=False)
                plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'_TR_IR', bbox_inches='tight')
                plt.show()
                
    ######################## Find maxima and plot peak wavenumber vs time #########
                singlet_peaks = []
                triplet_peaks = []
                
                for i, td in enumerate(time_delays):
                    j = 0
                    for wn in heat_map.iloc[td][20:107]:

                        if wn > max(heat_map.iloc[td][20:107])/3:
                            if td < 50:
                                singlet_peaks.append(float(j+20))
                            else:
                                triplet_peaks.append(float(j+20))
                        j += 1
                    
                    j = 0
                    for wn in heat_map.iloc[td][128:255]:
                        
                        if wn > 2*max(heat_map.iloc[td][128:255])/3:
                            if td < 50:
                                singlet_peaks.append(float(j+128))
                            else:
                                triplet_peaks.append(float(j+128))
                        j += 1
                    
                singlet_peaks = sorted(singlet_peaks)
                singlet_peakss = []
                singlet_peaksss = []
                for sp in singlet_peaks:
                    if sp not in singlet_peakss:
                        singlet_peakss.append(sp)

                for sp in singlet_peakss:
                    if math.isnan(np.average(heat_map.iloc[:][sp])):
                        continue
                    else:
                        singlet_peaksss.append(sp)

                singlet_peaks = [int(sps) for sps in singlet_peaksss]

                singlet_peaks = singlet_peaks[::3]
                #singlet_peaks = singlet_peaks[1:]

                triplet_peaks = sorted(triplet_peaks)
                triplet_peakss = []
                for sp in triplet_peaks:
                    if sp not in singlet_peakss:
                        triplet_peakss.append(sp)
                triplet_peaks = [int(sps) for sps in triplet_peakss]
                
                
################################################################################                
                
                fig,(ax,ax2) = plt.subplots(1, 2, sharey=False)            
                # plot the same data on both axes
                
                for td in time_delays:
                    ax.plot(wavenumbers[20:107], heat_map.iloc[td][20:107])
                    ax2.plot(wavenumbers[128:255], heat_map.iloc[td][128:255], label='time='+str(time[td])+'ps')
                ax.vlines(x=wavenumbers[[k for k in singlet_peaks if k <= 107]], ymin = min(heat_map.iloc[time_delays[2]][20:107]), ymax = max(heat_map.iloc[time_delays[2]][20:107]), color='black', linestyle='--')
                ax2.vlines(x=wavenumbers[[k for k in singlet_peaks if k >= 128]], ymin = min(heat_map.iloc[time_delays[2]][128:255]), ymax = max(heat_map.iloc[time_delays[2]][128:255]), color='black', linestyle='--')
                ax.yaxis.tick_left()
                ax2.yaxis.tick_right()
                ax.set_xlabel('Wavenumbers / cm$^{-1}$')
                ax2.set_xlabel('Wavenumbers / cm$^{-1}$')
                ax.set_ylabel(r'$\Delta$ Absorbance')
                plt.title(filename)
                plt.legend(loc='center right', bbox_to_anchor=(2.2, 0.5),
                           ncol=1, fancybox=True, shadow=False)
                #plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'_TR_IR', bbox_inches='tight')
                plt.show()
                
###############################################################################                
                
                # Curve fit to extract lifetime
                x = []
                y = []
                for sp in singlet_peaks:
                    filenames_4_df.append(filename)
                    peak_4_df_pixel.append(sp)
                    peak_4_df.append(wavenumbers[sp])
                    max_peak_time = np.where(heat_map.iloc[:] == max(heat_map.iloc[:64][sp]))[0][0]
                    xdata = time[max_peak_time:].astype(float)
                    xdata = xdata[:52]
                    x.append(xdata - xdata[0])
                    y.append(heat_map.iloc[:][sp].to_numpy().transpose()[max_peak_time:max_peak_time+52])
                     

                # #MONO, SYMFIT
                # B = sf.Parameter('B',value=0.0005, min=0.0, max=1.0)
                # iy = 0
                # As = []
                # Cs = []
                # xs = []
                # ys = []
                # for i, ffs in enumerate(singlet_peaks):
                #     iy = iy + 1
                #     globals()['y_%s' % iy], globals()['x_%s' % iy] = sf.variables(f'y_{iy},  x_{iy}')
                #     ys.append(globals()['y_%s' % iy])
                #     xs.append(globals()['x_%s' % iy])
                #     globals()['A_%s' % iy], globals()['C_%s' % iy] = sf.parameters(f'A_{iy}, C_{iy}')
                #     As.append(globals()['A_%s' % iy])
                #     Cs.append(globals()['C_%s' % iy])

                # # model_mono = []
                # # #model_mono = sf.Model({y_1: C_1 + A_1 * exp( - B * x_1)})
                # # iy = 0
                # # for i, ffs in enumerate(singlet_peaks):
                # #     iy = iy + 1
                # #     #model_mono.append(sf.Model({globals()['y_%s: C_%s + A_%s * exp(- B * x_%s' % {iy, iy, iy, iy}]}))
                # #     #model_mono.append({f'y_{iy+1}': f'C_{iy+1}' + f'A_{iy+1}' * exp(- B * f'x_{iy+1}')})
                # #     model_mono.append({ys[i]: Cs[i] + As[i] * exp(- B * xs[i])})
                # # model_monodd = sf.Model(model_mono)
                
                # model_mono = sf.Model({ys[i]: Cs[i] + As[i] * exp(- B * xs[i]) for i, fss in enumerate(ys)})
                
                # #model_mono = sf.Model({y_1: C_1 + A_1 * exp( - B * x_1),y_2: C_2 + A_2 * exp( - B * x_2)})

                # #fit = Fit(model=model_mono, x[i], y[i] for i, fss in enumerate(ys))
                # fit = Fit(model_mono, *x, *y)
                # fit_result = fit.execute()
                # fit_params_mono = fit_result.params
                # fit_error_mono = fit_result.covariance_matrix
                # fit_error_mono = np.sqrt(np.diag(fit_error_mono))
                # print(llllll)
                # print(fit_params_mono['A_1'])
                # print(fit_params_mono['B'])
                # print(fit_params_mono['C_1'])
                

                # plt.figure()
                # i = 0
                # for A in singlet_peaks:
                #     i = i + 1
                #     #plt.plot(x[i-1], Lifetime(x[i-1], fit_params_mono[f'A_{i}'], fit_params_mono['B'], fit_params_mono[f'C_{i}']))
                #     plt.plot(x[i-1], y[i-1], label=str(round(wavenumbers[A],1))+'cm$^{-1}$')
                #     plt.yscale('log')
                # plt.legend(loc='center right', bbox_to_anchor=(1.4, 0.5),
                #            ncol=1, fancybox=True, shadow=False)
                # plt.show()
                
                x = np.array(x)
                y = np.array(y)
                
                #MONO LMFIT
                fit_params_mono = Parameters()
                fit_params_mono.add('B', value=0.002, min=0.00001, max=0.5)
                for iy, ffs in enumerate(singlet_peaks):
                    fit_params_mono.add(f'A_{iy+1}', value=max(heat_map.iloc[:][ffs]), min=0.0, max=0.01)
                    fit_params_mono.add(f'C_{iy+1}', value=0.0, min=-0.01, max=0.01)

                out_mono = minimize(objectiveMono, fit_params_mono, args=(x, y))
                parameters_mono = out_mono.params
                errors_mono = out_mono.uvars
                
                # plt.figure()
                # i = 0
                # for A in singlet_peaks:
                #     i = i + 1
                #     plt.plot(x[i-1], Lifetime(x[i-1], parameters_mono[f'A_{i}'], parameters_mono['B'], parameters_mono[f'C_{i}']))
                #     plt.plot(x[i-1], y[i-1], label=str(round(wavenumbers[A],1))+'cm$^{-1}$')
                #     #plt.yscale('log')
                #     plt.xlim(-0.1, 10000)
                # plt.legend(loc='center right', bbox_to_anchor=(1.4, 0.5),
                #            ncol=1, fancybox=True, shadow=False)
                # plt.show()
                #Mono
                for i, sp in enumerate(singlet_peaks):
                    fit_A_mono.append(parameters_mono[f'A_{i+1}'].value)
                    fit_B_mono.append(parameters_mono['B'].value)
                    fit_C_mono.append(parameters_mono[f'C_{i+1}'].value)
                    
                    fit_A_mono_err.append(parameters_mono[f'A_{i+1}'].stderr)
                    fit_B_mono_err.append(parameters_mono['B'].stderr)
                    fit_C_mono_err.append(parameters_mono[f'C_{i+1}'].stderr)
                
                    chi_sq_mono.append(out_mono.redchi)
                    lifetime_mono.append(1/parameters_mono['B'].value)
                    try: 
                        lifetime_mono_err.append(np.sqrt((1/parameters_mono['B'].value)**2 - (1/(parameters_mono['B'].value + parameters_mono['B']))**2))
                    except ValueError:
                        lifetime_mono_err.append(np.nan)
                        continue
                    except TypeError:
                        lifetime_mono_err.append(np.nan)
                        continue 
                #DOUBLE
                
                fit_params_bi = Parameters()
                fit_params_bi.add('B', value=0.0005, min=0.00001, max=0.1)
                fit_params_bi.add('D', value=0.002, min=0.00001, max=0.5)
                for iy, fss in enumerate(singlet_peaks):
                    fit_params_bi.add(f'A_{iy+1}', value=max(heat_map.iloc[:][sp])/2, min=0.0, max=0.01)
                    fit_params_bi.add(f'C_{iy+1}', value=max(heat_map.iloc[:][sp])/2, min=0.0, max=0.01)
                    fit_params_bi.add(f'E_{iy+1}', value=0.0, min=-0.01, max=0.01)
                    
                out_bi = minimize(objectiveDouble, fit_params_bi, args=(x, y))
                parameters = out_bi.params
                
                #Double
                for i, sp in enumerate(singlet_peaks):
                    fit_A.append(parameters[f'A_{i+1}'].value)
                    fit_B.append(parameters['B'].value)
                    fit_C.append(parameters[f'C_{i+1}'].value)
                    fit_D.append(parameters['D'].value)
                    fit_E.append(parameters[f'E_{i+1}'].value)
                    
                    fit_A_err.append(parameters[f'A_{i+1}'].stderr)
                    fit_B_err.append(parameters['B'].stderr)
                    fit_C_err.append(parameters[f'C_{i+1}'].stderr)
                    fit_D_err.append(parameters['D'].stderr)
                    fit_E_err.append(parameters[f'E_{i+1}'].stderr)
                
                    lifetime1.append(1/parameters['B'].value)
                    lifetime2.append(1/parameters['D'].value)
                    av_amp_lifetime.append(((parameters[f'A_{i+1}'].value*(1/parameters['B'].value))+(parameters[f'C_{i+1}'].value*(1/parameters['D'].value)))/(parameters[f'A_{i+1}'].value+parameters[f'C_{i+1}'].value))
                    av_int_lifetime.append(((parameters[f'A_{i+1}'].value*((1/parameters['B'].value)**2))+(parameters[f'C_{i+1}'].value*((1/parameters['D'].value)**2)))/((parameters[f'A_{i+1}'].value*(1/parameters['B'].value))+(parameters[f'C_{i+1}'].value*(1/parameters['D'].value))))
                    
                    chi_sq.append(out_bi.chisqr)

                    try:
                        lifetime1_err.append(np.sqrt((1/parameters['B'].value)**2 - (1/(parameters['B'].value + parameters['B'].stderr))**2))
                        lifetime2_err.append(np.sqrt((1/parameters['D'].value)**2 - (1/(parameters['D'].value + parameters['D'].stderr))**2))
                    except ValueError:
                        lifetime1_err.append(np.nan)
                        lifetime2_err.append(np.nan)
                        continue
                    except TypeError:
                        lifetime1_err.append(np.nan)
                        lifetime2_err.append(np.nan)
                        continue
                #TRIPLE
                fit_params_trip = Parameters()
                fit_params_trip.add('B', value=0.0005, min=0.00001, max=0.03)
                fit_params_trip.add('D', value=0.002, min=0.00001, max=0.1)
                fit_params_trip.add('F', value=0.005, min=0.00001, max=0.5)
                for iy, ffs in enumerate(singlet_peaks):
                    fit_params_trip.add(f'A_{iy+1}', value=max(heat_map.iloc[:][sp])/3, min=0.0, max=0.01)
                    fit_params_trip.add(f'C_{iy+1}', value=max(heat_map.iloc[:][sp])/3, min=0.0, max=0.01)
                    fit_params_trip.add(f'E_{iy+1}', value=max(heat_map.iloc[:][sp])/3, min=0.0, max=0.01)
                    fit_params_trip.add(f'G_{iy+1}', value=0.0, min=-0.01, max=0.01)
                    
                out_trip = minimize(objectiveTriple, fit_params_trip, args=(x, y))
                parameters_trip = out_trip.params
                
                #Triple
                for i, sp in enumerate(singlet_peaks):
                    fit_A_trip.append(parameters_trip[f'A_{i+1}'].value)
                    fit_B_trip.append(parameters_trip['B'].value)
                    fit_C_trip.append(parameters_trip[f'C_{i+1}'].value)
                    fit_D_trip.append(parameters_trip['D'].value)
                    fit_E_trip.append(parameters_trip[f'E_{i+1}'].value)
                    fit_F_trip.append(parameters_trip['F'].value)
                    fit_G_trip.append(parameters_trip[f'G_{i+1}'].value)
                    
                    fit_A_trip_err.append(parameters_trip[f'A_{i+1}'].stderr)
                    fit_B_trip_err.append(parameters_trip['B'].stderr)
                    fit_C_trip_err.append(parameters_trip[f'C_{i+1}'].stderr)
                    fit_D_trip_err.append(parameters_trip['D'].stderr)
                    fit_E_trip_err.append(parameters_trip[f'E_{i+1}'].stderr)
                    fit_F_trip_err.append(parameters_trip['F'].stderr)
                    fit_G_trip_err.append(parameters_trip[f'G_{i+1}'].stderr)
                    
                    lifetime1_trip.append(1/parameters_trip['B'].value)
                    lifetime2_trip.append(1/parameters_trip['D'].value)
                    lifetime3_trip.append(1/parameters_trip['F'].value)
                    
                    av_amp_lifetime_trip.append(((parameters_trip[f'A_{i+1}'].value*(1/parameters_trip['B'].value))+(parameters_trip[f'C_{i+1}'].value*(1/parameters_trip['D'].value))+(parameters_trip[f'E_{i+1}'].value*(1/parameters_trip['F'].value)))/(parameters_trip[f'A_{i+1}'].value+parameters_trip[f'C_{i+1}'].value+parameters_trip[f'E_{i+1}'].value))
                    av_int_lifetime_trip.append(((parameters_trip[f'A_{i+1}'].value*((1/parameters_trip['B'].value)**2))+(parameters_trip[f'C_{i+1}'].value*((1/parameters_trip['D'].value)**2))+(parameters_trip[f'E_{i+1}'].value*((1/parameters_trip['F'].value)**2)))/((parameters_trip[f'A_{i+1}'].value*(1/parameters_trip['B'].value))+(parameters_trip[f'C_{i+1}'].value*(1/parameters_trip['D'].value))+(parameters_trip[f'E_{i+1}'].value*(1/parameters_trip['F'].value))))
                
                    chi_sq_trip.append(out_trip.chisqr)
                    try:
                        lifetime1_trip_err.append(np.sqrt((1/parameters_trip['B'].value)**2 - (1/(parameters_trip['B'].value + parameters_trip['B'].stderr))**2))
                        lifetime2_trip_err.append(np.sqrt((1/parameters_trip['D'].value)**2 - (1/(parameters_trip['D'].value + parameters_trip['D'].stderr))**2))
                        lifetime3_trip_err.append(np.sqrt((1/parameters_trip['F'].value)**2 - (1/(parameters_trip['F'].value + parameters_trip['F'].stderr))**2))
                    
                    except ValueError:
                        lifetime1_trip_err.append(np.nan)
                        lifetime2_trip_err.append(np.nan)
                        lifetime3_trip_err.append(np.nan)
                        continue
                    except TypeError:
                        lifetime1_trip_err.append(np.nan)
                        lifetime2_trip_err.append(np.nan)
                        lifetime3_trip_err.append(np.nan)
                        continue
                
                i = 0
                for A in singlet_peaks:
                    fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, height_ratios={0.75,0.25})
                    i = i + 1
                    y_pred_mono = Lifetime(x[i-1], parameters_mono[f'A_{i}'], parameters_mono['B'], parameters_mono[f'C_{i}'])
                    y_pred_bi = LifetimeDoubleExp(x[i-1], parameters[f'A_{i}'], parameters['B'], 
                                                       parameters[f'C_{i}'], parameters['D'], parameters[f'E_{i}'])
                    y_pred_tri = LifetimeTripleExp(x[i-1], parameters_trip[f'A_{i}'], parameters_trip['B'], 
                                                       parameters_trip[f'C_{i}'], parameters_trip['D'], parameters_trip[f'E_{i}'],
                                                       parameters_trip['F'], parameters_trip[f'G_{i}'])
                    
                    axs[0].plot(x[i-1], y[i-1], label=str(round(wavenumbers[A],1))+'cm$^{-1}$')
                    axs[0].plot(x[i-1], y_pred_mono, linestyle='--')
                    axs[0].plot(x[i-1], y_pred_bi, linestyle='--')
                    axs[0].plot(x[i-1], y_pred_tri, linestyle='--')
                    axs[1].set_xlabel('Time Delay / ps')
                    axs[0].set_ylabel(r'$\Delta$ Absorbance')
                    plt.xlim(-0.1, 10000)
                    axs[0].set_title(compound+' '+solvent+' '+str(round(wavenumbers[A],1))+'cm-1')
                    
                    axs[1].plot((0,0))
                    axs[1].set_ylabel('Residuals')
                    axs[1].plot(x[i-1], y[i-1]-y_pred_mono, label='Mono')
                    axs[1].plot(x[i-1], y[i-1]-y_pred_bi, label='Double')
                    axs[1].plot(x[i-1], y[i-1]-y_pred_tri, label='Triple')
                    #axs[1].hlines(0, xdata[0]-(xdata[-1]*0.06), xdata[-1]*1.06, colors='Black', linestyles='solid',alpha=0.4)
                    #fig.set_major_locator(plt.MaxNLocator(10))
                    #axs[1].set_xlim(xdata[0]-(xdata[-1]*0.05), xdata[-1]*1.05)
                    plt.locator_params(axis='x',nbins=10)
                    axs[0].legend()
                    print(x[i-1][-1])
                    if x[i-1][-1] > 10000:
                        plt.xlim(-200, 10000)
                        axs[1].hlines(0, -200, 10000, colors='Black', linestyles='solid',alpha=0.4)
                    else:
                        plt.xlim(x[i-1][0]-(x[i-1][-1]*0.05), x[i-1][-1]*1.05)
                        axs[1].hlines(0, x[i-1][0]-(x[i-1][-1]*0.06), x[i-1][-1]*1.06, colors='Black', linestyles='solid',alpha=0.4)
                    plt.legend(loc='center right', bbox_to_anchor=(1.4, 0.5),
                               ncol=1, fancybox=True, shadow=False)
                    plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+' '+str(round(wavenumbers[A]))+'cm-1', bbox_inches='tight')
                    plt.show()
                    
                    residuals_mono.append(sum(abs(y[i-1]-y_pred_mono)))
                    residuals.append(sum(abs(y[i-1]-y_pred_bi)))
                    residuals_trip.append(sum(abs(y[i-1]-y_pred_tri)))

                # for sp in singlet_peaks[i]:
                #     try:

                #     except RuntimeError:
                #         continue
                #     except ValueError:
                #         continue
                #     except TypeError:
                #         continue
            
            # plt.plot(time, heat_map.iloc[:][max_peak], label='Measured @ '+str(round(wavenumbers[max_peak]))+'cm-1')
            # y_pred = Lifetime(time.to_numpy().astype(float), *parameters)
            # plt.plot(time, y_pred, label=r'Best Fit, $\tau$='+ str(round(1/parameters[1]))+'ps')
            # plt.xlabel('Time Delay / ps')
            # plt.ylabel(r'$\Delta$ Absorbance')
            # plt.title(filename)
            # #fig.set_major_locator(plt.MaxNLocator(10))
            # plt.locator_params(axis='x',nbins=10)
            # plt.legend()
            # #plt.savefig(Results_directory+'\\'+compound+'\\'+solvent+'\\'+str(filename)+'lifetime_difxscale')
            # plt.show()
                    


        df_lifetime = pd.DataFrame({'filename': filenames_4_df,
                                    'Peak position': peak_4_df_pixel, 'Peak position / cm-1': peak_4_df, 
                                    'fit_A_mono':fit_A_mono, 'fit_A_mono_err':fit_A_mono_err,
                                    'fit_B_mono':fit_B_mono, 'fit_B_mono_err':fit_B_mono_err,
                                    'fit_C_mono':fit_C_mono, 'fit_C_mono_err':fit_C_mono_err,
                                    'chi_sq_mono': chi_sq_mono, 'residuals_mono': residuals_mono,
                                    'lifetime mono':lifetime_mono, 'lifetime_mono_err':lifetime_mono_err,
                                    'fit_A':fit_A, 'fit_A_err':fit_A_err,
                                    'fit_B':fit_B, 'fit_B_err':fit_B_err,
                                    'fit_C':fit_C, 'fit_C_err':fit_C_err,
                                    'fit_D':fit_D, 'fit_D_err':fit_D_err,
                                    'fit_E':fit_E, 'fit_E_err':fit_E_err,
                                    'residuals_bi': residuals, 'chi_sq_bi': chi_sq,
                                    'Bi lifetime tau1':lifetime1, 'Bi lifetime_tau1_err':lifetime1_err,
                                    'Bi lifetime tau2':lifetime2, 'Bi lifetime_tau2_err':lifetime2_err,
                                    'Bi Amp Average lifetime / ps': av_amp_lifetime, 'Bi Intensity Average lifetime / ps': av_int_lifetime,
                                    'fit_A_trip':fit_A_trip, 'fit_A_trip_err':fit_A_trip_err,
                                    'fit_B_trip':fit_B_trip, 'fit_B_trip_err':fit_B_trip_err,
                                    'fit_C_trip':fit_C_trip, 'fit_C_trip_err':fit_C_trip_err,
                                    'fit_D_trip':fit_D_trip, 'fit_D_trip_err':fit_D_trip_err,
                                    'fit_E_trip':fit_E_trip, 'fit_E_trip_err':fit_E_trip_err,
                                    'fit_F_trip':fit_F_trip, 'fit_F_trip_err':fit_F_trip_err,
                                    'fit_G_trip':fit_G_trip, 'fit_G_trip_err':fit_G_trip_err,
                                    'chi_sq_trip': chi_sq_trip, 'residuals_trip': residuals_trip,
                                    'Trip lifetime tau1':lifetime1_trip, 'Trip lifetime_tau1_err':lifetime1_trip_err,
                                    'Trip lifetime tau2':lifetime2_trip, 'Trip lifetime_tau2_err':lifetime2_trip_err,
                                    'Trip lifetime tau3':lifetime3_trip, 'Trip lifetime_tau3_err':lifetime3_trip_err,
                                    'Trip Amp Average lifetime / ps': av_amp_lifetime_trip, 'Trip Intensity Average lifetime / ps': av_int_lifetime_trip})
        
        df_list.append(df_lifetime)
        
    Excelwriter = pd.ExcelWriter(Results_directory+'\\'+compound+' Singlet Lifetimes.xlsx',engine="xlsxwriter")
    for i, df in enumerate (df_list):    
        df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
        
    Excelwriter.save()

    
    
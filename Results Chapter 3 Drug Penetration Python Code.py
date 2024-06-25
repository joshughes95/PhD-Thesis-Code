# -*- coding: utf-8 -*-
"""
Created on Thu 01 Dec 22 13:43

@author: tscw79
"""

from os.path import dirname, join as pjoin
import scipy.io as sio
import numpy as np
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

############################ MANUAL IMPUTS ####################################
pixel_size = 0.32933005992010652463382157123835

input_path = askdirectory(title='Select Input Directory')
input_path_control = askdirectory(title='Select Control Directory')
output_path = askdirectory(title='Select Output Directory')

############################ Curve Fitting ####################################

def Uptake(x, A, B, C):
    y = A*np.exp(-B*x)+C
    return y

def Diffusion(x, M0, D, t, C):
    y = M0*(1-np.erf(x/(np.sqrt(4*D*t))))+C
    return y

def Alignment(x, M, C):
    y = M*x + C
    return y

############################ empty list to collate dataframes #################
repeats = len(input_path)
df_list = []
sheet_names = []

integrated_fluorescence_epi = []
integrated_fluorescence_err_epi = []

integrated_fluorescence_der = []
integrated_fluorescence_err_der = []

integrated_fluorescence_cf_epi = []
integrated_fluorescence_cf_err_epi = []

integrated_fluorescence_cf_der = []
integrated_fluorescence_cf_err_der = []

########################### CONTROLS #########################################


epidermis_sizes = []
df_control = pd.DataFrame()
    
for j, filename in enumerate(os.listdir(input_path_control)):
    
   img = io.imread(input_path_control+'\\'+filename) # open image
   control_img = img[0]   
   mask_img = img[1]                                       # select 129 channel  
        
   img_2_view = (control_img/256).astype('uint8')
   c = 255 / np.log(1 + np.max(img_2_view))
   log_image = (c * (np.log(img_2_view + 1))).astype('uint8')
   plt.imshow(log_image)                                 # view image
   plt.colorbar()
   plt.show()
    
############################# Adjust constrast and view image #################
        # plt.hist(img.flat, bins=100, range=(0,255))
        # equ = cv2.equalizeHist(img)
        # #see equalised histogram
        # plt.hist(equ.flat, bins=100, range=(0,255))
        # plt.show()
        # plt.imshow(equ)
        # plt.show()
        
############################## Remove 0's to align top of image ###############
   control_epidermis = []
   control_dermis = []
   ydata = []      # array to hold length of zeros before start of tissue
   epidermis_size = []
   control_transpose = np.transpose(control_img)
   mask_transpose = np.transpose(mask_img)
   # img2 = img2[10:len(img2)-10] # removes stray pixels not picked up by the magic wand tool

   for i, (control_col, mask_col) in enumerate(zip(control_transpose, mask_transpose)):
        control_epidermis.append(control_col[mask_col==1])
        control_dermis.append(control_col[mask_col==2])
        ydata.append(len(control_col[mask_col==0]))
        epidermis_size.append(len(control_col[mask_col==1]))
        
   epidermis_size = np.mean(epidermis_size)
        
   # calculate slope of misalignment
   xdata = np.linspace(0, len(ydata)*pixel_size, len(ydata)).reshape((-1, 1))
   model = LinearRegression().fit(xdata, ydata)
   gradient = model.coef_
   intercept = model.intercept_
   r_sq = model.score(xdata, ydata)
        
   #epidermis
   transposed_tuples_epi = list(zip(*control_epidermis))
   transposed_epi = [list(sublist) for sublist in transposed_tuples_epi]
   lll_epi = transposed_epi
        
   tranposed_tuples_epi = zip_longest(*control_epidermis, fillvalue=None)
   transposed_tuples_list_epi = list(tranposed_tuples_epi)
   transposed_epi = [list(sublist) for sublist in transposed_tuples_list_epi]
   transposed_epi = np.array(transposed_epi)
    
   #dermis
   transposed_tuples_der = list(zip(*control_dermis))
   transposed_der = [list(sublist) for sublist in transposed_tuples_der]
   lll_der = transposed_der
    
   tranposed_tuples_der = zip_longest(*control_dermis, fillvalue=None)
   transposed_tuples_list_der = list(tranposed_tuples_der)
   transposed_der = [list(sublist) for sublist in transposed_tuples_list_der]
   transposed_der = np.array(transposed_der)
    
############################### Calculate mean & std along x-axis #############    
   #epidermis
   mean_epi = np.zeros(len(transposed_epi))
   std_epi = np.zeros(len(transposed_epi))
   for i, col in enumerate(transposed_epi):
       res_epi = list(filter(lambda item: item is not None, col))
       mean_epi[i] = np.mean(res_epi)
       std_epi[i] = np.std(res_epi)
            
   #dermis
   mean_der = np.zeros(len(transposed_der))
   std_der = np.zeros(len(transposed_der))
   for i, col in enumerate(transposed_der):
       res_der = list(filter(lambda item: item is not None, col))
       mean_der[i] = np.mean(res_der)
       std_der[i] = np.std(res_der)
        
############################### Add to Dataframe ##############################
   conversion_factor = np.cos(np.arctan(gradient))
   epidermis_sizes.append(epidermis_size*conversion_factor)
           ############### CHECK IF THIS IS TRUE PIXEL SIZE #########################
   x_axis_epi = np.arange(0, len(mean_epi), 1)*conversion_factor*pixel_size
   x_axis_der = np.arange(0, len(mean_der), 1)*conversion_factor*pixel_size + x_axis_epi[-1]
   additional = pd.DataFrame({'Depth Repeat'+str(j+1)+' Epi':x_axis_epi, 'Repeat'+str(j+1)+' Epi':mean_epi, 'Repeat_'+str(j+1)+' Epi_Std':std_epi})
   new = pd.concat([df_control, additional], axis=1)
   df_control = new
        
   additional = pd.DataFrame({'Depth Repeat'+str(j+1)+' Der':x_axis_der, 'Repeat'+str(j+1)+' Der':mean_der, 'Repeat_'+str(j+1)+' Der_Std':std_der})
   new = pd.concat([df_control, additional], axis=1)
   df_control = new
    
        
   #x_axis = np.arange(0, len(transposed), 1)*pixel_size
   plt.plot(x_axis_epi, mean_epi, label='Epidermis')
   plt.plot(x_axis_der, mean_der, label='Dermis')
   plt.title(filename+'129 depth penetration')
   plt.show()

############################### Average over Individual images ################
array_lengths_epi = []
x_values_epi = []
y_values_epi = []
    
array_lengths_der = []
x_values_der = []
y_values_der = []
for k in range(j+1):
    array_epi = df_control.iloc[:,6*k+1]
    array_epi = array_epi[~np.isnan(df_control.iloc[:,6*k+1])]
    array_epi = array_epi[array_epi>0]
    array_lengths_epi.append(len(array_epi))
    x_values_epi.append(df_control.iloc[:,6*k])
    y_values_epi.append(df_control.iloc[:,6*k+1])
        
    array_der = df_control.iloc[:,6*k+4]
    array_der = array_der[~np.isnan(df_control.iloc[:,6*k+4])]
    array_der = array_der[array_der>0]
    array_lengths_der.append(len(array_der))
    x_values_der.append(df_control.iloc[:,6*k+3] - df_control.iloc[0,6*k+3])
    y_values_der.append(df_control.iloc[:,6*k+4])
        
# interpolate arrays to same x_values
min_length_epi = min(array_lengths_epi)
min_length_der = min(array_lengths_der)
    
y_values_interpolated_epi = []
y_values_interpolated_epi.append(y_values_epi[0][:min_length_epi])
    
y_values_interpolated_der = []
y_values_interpolated_der.append(y_values_der[0][:min_length_der])
    
for k in range(len(x_values_epi)-1):
    y_values_interpolated_epi.append(np.interp(x_values_epi[0][:min_length_epi], 
                                               x_values_epi[k+1][:min_length_epi], 
                                               y_values_epi[k+1][:min_length_epi]))
        
for k in range(len(x_values_der)-1):
    y_values_interpolated_der.append(np.interp(x_values_der[0][:min_length_der], 
                                               x_values_der[k+1][:min_length_der], 
                                               y_values_der[k+1][:min_length_der]))
    
# Calculate mean and Std
Mean_epi = np.mean(y_values_interpolated_epi, axis=0)
Std_epi = np.std(y_values_interpolated_epi, axis=0)
    
############### CHECK IF THIS IS TRUE PIXEL SIZE #########################
x_axis_epi = np.arange(0, min_length_epi, 1)*x_values_epi[0][1]
    
Mean_der = np.mean(y_values_interpolated_der, axis=0)
Std_der = np.std(y_values_interpolated_der, axis=0)
    
epidermis_size_mean = np.mean(epidermis_sizes)
x_axis_der = np.arange(0, min_length_der, 1)*x_values_epi[0][1] + x_axis_epi[-1]
    
additional = pd.DataFrame({'Depth Epi': x_axis_epi, 'Mean Epi':Mean_epi, 'Std Epi':Std_epi})
new = pd.concat([df_control, additional], axis=1)
df_control = new
    
additional = pd.DataFrame({'Depth Der': x_axis_der, 'Mean Der':Mean_der, 'Std Der':Std_der})
new = pd.concat([df_control, additional], axis=1)
df_control = new

############################ SAVE CONTROL DATA ################################

Excelwriter = pd.ExcelWriter(output_path+'\\Control Results.xlsx',engine="xlsxwriter")
 
df_control.to_excel(Excelwriter, index=False)
    
Excelwriter.save()

############################ Control averages #################################

epidermis_control = np.mean(Mean_epi)
dermis_control = np.mean(Mean_der)
control_x_epi = x_axis_epi
control_y_epi = Mean_epi
control_x_der = x_axis_der
control_y_der = Mean_der

plt.plot(x_axis_epi, Mean_epi, label='Epidermis')
plt.plot(x_axis_der, Mean_der, label='Dermis')
plt.legend()
plt.title('Control')
plt.savefig(output_path+'//LightOx129 Control in Pig Tongue.png')
plt.show()

########################### REPEATS ###########################################

############################ loop over folders and individual images ##########

for k, folder in enumerate(os.listdir(input_path)):
    df = pd.DataFrame()
    sheet_names.append(folder)
    epidermis_sizes = []
    
    for j, filename in enumerate(os.listdir(input_path+'\\'+folder)):
    
        img = io.imread(input_path+'\\'+folder+'\\'+filename) # open image
        compound_img = img[0]   
        mask_img = img[1]                                       # select 129 channel  
        
        # img_2_view = (compound_img/256).astype('uint8')
        # c = 255 / np.log(1 + np.max(img_2_view))
        # log_image = (c * (np.log(img_2_view + 1))).astype('uint8')
        # plt.imshow(log_image)                                 # view image
        # plt.colorbar()
        # plt.show()
    
############################# Adjust constrast and view image #################
        # plt.hist(img.flat, bins=100, range=(0,255))
        # equ = cv2.equalizeHist(img)
        # #see equalised histogram
        # plt.hist(equ.flat, bins=100, range=(0,255))
        # plt.show()
        # plt.imshow(equ)
        # plt.show()
        
############################## Remove 0's to align top of image ###############
        compound_epidermis = []
        compound_dermis = []
        ydata = []      # array to hold length of zeros before start of tissue
        epidermis_size = []
        compound_transpose = np.transpose(compound_img)
        mask_transpose = np.transpose(mask_img)
        # img2 = img2[10:len(img2)-10] # removes stray pixels not picked up by the magic wand tool

        for i, (compound_col, mask_col) in enumerate(zip(compound_transpose, mask_transpose)):
            compound_epidermis.append(compound_col[mask_col==1])
            compound_dermis.append(compound_col[mask_col==2])
            ydata.append(len(compound_col[mask_col==0]))
            epidermis_size.append(len(compound_col[mask_col==1]))
        
        epidermis_size = np.mean(epidermis_size)
        
        # calculate slope of misalignment
        xdata = np.linspace(0, len(ydata)*pixel_size, len(ydata)).reshape((-1, 1))
        model = LinearRegression().fit(xdata, ydata)
        gradient = model.coef_
        intercept = model.intercept_
        r_sq = model.score(xdata, ydata)
        
        #epidermis
        transposed_tuples_epi = list(zip(*compound_epidermis))
        transposed_epi = [list(sublist) for sublist in transposed_tuples_epi]
        lll_epi = transposed_epi
        
        tranposed_tuples_epi = zip_longest(*compound_epidermis, fillvalue=None)
        transposed_tuples_list_epi = list(tranposed_tuples_epi)
        transposed_epi = [list(sublist) for sublist in transposed_tuples_list_epi]
        transposed_epi = np.array(transposed_epi)
        
        #dermis
        transposed_tuples_der = list(zip(*compound_dermis))
        transposed_der = [list(sublist) for sublist in transposed_tuples_der]
        lll_der = transposed_der
        
        tranposed_tuples_der = zip_longest(*compound_dermis, fillvalue=None)
        transposed_tuples_list_der = list(tranposed_tuples_der)
        transposed_der = [list(sublist) for sublist in transposed_tuples_list_der]
        transposed_der = np.array(transposed_der)
    
############################### Calculate mean & std along x-axis #############    
        #epidermis
        mean_epi = np.zeros(len(transposed_epi))
        std_epi = np.zeros(len(transposed_epi))
        for i, col in enumerate(transposed_epi):
            res_epi = list(filter(lambda item: item is not None, col))
            mean_epi[i] = np.mean(res_epi)
            std_epi[i] = np.std(res_epi)
            
        #dermis
        mean_der = np.zeros(len(transposed_der))
        std_der = np.zeros(len(transposed_der))
        for i, col in enumerate(transposed_der):
            res_der = list(filter(lambda item: item is not None, col))
            mean_der[i] = np.mean(res_der)
            std_der[i] = np.std(res_der)
        
############################### Add to Dataframe ##############################
        conversion_factor = np.cos(np.arctan(gradient))
        epidermis_sizes.append(epidermis_size*conversion_factor)
            ############### CHECK IF THIS IS TRUE PIXEL SIZE #########################
        x_axis_epi = np.arange(0, len(mean_epi), 1)*conversion_factor*pixel_size
        x_axis_der = np.arange(0, len(mean_der), 1)*conversion_factor*pixel_size + x_axis_epi[-1]
        additional = pd.DataFrame({'Depth Repeat'+str(j+1)+' Epi':x_axis_epi, 'Repeat'+str(j+1)+' Epi':mean_epi, 'Repeat_'+str(j+1)+' Epi_Std':std_epi})
        new = pd.concat([df, additional], axis=1)
        df = new
        
        additional = pd.DataFrame({'Depth Repeat'+str(j+1)+' Der':x_axis_der, 'Repeat'+str(j+1)+' Der':mean_der, 'Repeat_'+str(j+1)+' Der_Std':std_der})
        new = pd.concat([df, additional], axis=1)
        df = new
    
        
        #x_axis = np.arange(0, len(transposed), 1)*pixel_size
        plt.plot(x_axis_epi, mean_epi, label='Epidermis')
        plt.plot(x_axis_der, mean_der, label='Dermis')
        plt.title(filename+'129 depth penetration')
        plt.show()

############################### Average over Individual images ################
    array_lengths_epi = []
    x_values_epi = []
    y_values_epi = []
    
    array_lengths_der = []
    x_values_der = []
    y_values_der = []
    for i in range(j+1):
        array_epi = df.iloc[:,6*k+1]
        array_epi = array_epi[~np.isnan(df.iloc[:,6*i+1])]
        array_epi = array_epi[array_epi>0]
        array_lengths_epi.append(len(array_epi))
        x_values_epi.append(df.iloc[:,6*i])
        y_values_epi.append(df.iloc[:,6*i+1])
        
        array_der = df.iloc[:,6*k+4]
        array_der = array_der[~np.isnan(df.iloc[:,6*i+4])]
        array_der = array_der[array_der>0]
        array_lengths_der.append(len(array_der))
        x_values_der.append(df.iloc[:,6*i+3] - df.iloc[0,6*i+3])
        y_values_der.append(df.iloc[:,6*i+4])
        
    # interpolate arrays to same x_values
    min_length_epi = min(array_lengths_epi)
    min_length_der = min(array_lengths_der)
    
    y_values_interpolated_epi = []
    y_values_interpolated_epi_control = np.interp(x_values_epi[0][:min_length_epi], 
                                                  control_x_epi[:min_length_epi], 
                                                  control_y_epi[:min_length_epi])
        
    y_values_interpolated_epi.append(y_values_epi[0][:min_length_epi])
    
    y_values_interpolated_der = []
    y_values_interpolated_der_control = np.interp(x_values_der[0][:min_length_der], 
                                                  control_x_der[:min_length_der], 
                                                  control_y_der[:min_length_der])
    y_values_interpolated_der.append(y_values_der[0][:min_length_der])
    
    for i in range(len(x_values_epi)-1):
        y_values_interpolated_epi.append(np.interp(x_values_epi[0][:min_length_epi], 
                                                   x_values_epi[i+1][:min_length_epi], 
                                                   y_values_epi[i+1][:min_length_epi]) - y_values_interpolated_epi_control)
        
    for i in range(len(x_values_der)-1):
        y_values_interpolated_der.append(np.interp(x_values_der[0][:min_length_der], 
                                                   x_values_der[i+1][:min_length_der], 
                                                   y_values_der[i+1][:min_length_der]) - y_values_interpolated_der_control)
    
    # Calculate mean and Std
    Mean_epi = np.mean(y_values_interpolated_epi, axis=0)
    Std_epi = np.std(y_values_interpolated_epi, axis=0)
    
    ############### CHECK IF THIS IS TRUE PIXEL SIZE #########################
    x_axis_epi = np.arange(0, min_length_epi, 1)*x_values_epi[0][1]
    
    Mean_der = np.mean(y_values_interpolated_der, axis=0)
    Std_der = np.std(y_values_interpolated_der, axis=0)
    
    epidermis_size_mean = np.mean(epidermis_sizes)
    x_axis_der = np.arange(0, min_length_der, 1)*x_values_epi[0][1] + x_axis_epi[-1]
    
    additional = pd.DataFrame({'Depth Epi': x_axis_epi, 'Mean Epi':Mean_epi, 'Std Epi':Std_epi})
    new = pd.concat([df, additional], axis=1)
    df = new
    
    additional = pd.DataFrame({'Depth Der': x_axis_der, 'Mean Der':Mean_der, 'Std Der':Std_der})
    new = pd.concat([df, additional], axis=1)
    df = new
    
################################## Curve Fitting #############################
    #Epidermis
    xdata_epi = x_axis_epi[10:] 

    ydata_epi = Mean_epi[10:]
    ydata_err_epi = Std_epi[10:]
    
    parameters_epi, covariance_epi = curve_fit(Uptake, xdata_epi, ydata_epi)
    y_pred_epi = Uptake(xdata_epi, *parameters_epi)
    r_sq_epi = r2_score(ydata_epi, y_pred_epi)
    errors_epi = np.sqrt(np.diag(covariance_epi))

    fit_A_epi = parameters_epi[0]
    fit_B_epi = parameters_epi[1]
    fit_C_epi = parameters_epi[2]
    
    fit_A_err_epi = errors_epi[0]
    fit_B_err_epi = errors_epi[1]
    fit_C_err_epi = errors_epi[2]
    
    #Dermis
    xdata_der = x_axis_der[:-500] - x_axis_der[0]

    ydata_der = Mean_der[:-500]
    ydata_err_der = Std_der
    
    plt.plot(xdata_der, ydata_der)
    plt.show()
    
    parameters_der, covariance_der = curve_fit(Uptake, xdata_der, ydata_der, p0=[250,1,250], maxfev=5000, bounds=((0,0,0),(np.inf,np.inf,np.inf)))
    y_pred_der = Uptake(xdata_der, *parameters_der)
    r_sq_der = r2_score(ydata_der, y_pred_der)
    errors_der = np.sqrt(np.diag(covariance_der))

    fit_A_der = parameters_der[0]
    fit_B_der = parameters_der[1]
    fit_C_der = parameters_der[2]
    
    fit_A_err_der = errors_der[0]
    fit_B_err_der = errors_der[1]
    fit_C_err_der = errors_der[2]
    
    additional = pd.DataFrame({'fit_A_Epi':[fit_A_epi], 'fit_A_err_Epi':[fit_A_err_epi],
                               'fit_B_Epi':[fit_B_epi], 'fit_B_err_Epi':[fit_B_err_epi],
                               'fit_C_Epi':[fit_C_epi], 'fit_C_err_Epi':[fit_C_err_epi],'r_sq_Epi': [r_sq_epi],
                               'fit_A_Der':[fit_A_der], 'fit_A_err_Der':[fit_A_err_der],
                               'fit_B_Der':[fit_B_der], 'fit_B_err_Der':[fit_B_err_der],
                               'fit_C_Der':[fit_C_der], 'fit_C_err_Der':[fit_C_err_der],'r_sq_Der': [r_sq_der]})
    
    new = pd.concat([df, additional], axis=1)
    df = new
    df_list.append(df)
    
    int_fluo_sums_epi = np.sum(y_values_interpolated_epi, axis=1)*x_axis_epi[1]
    int_fluo_sums_der = np.sum(y_values_interpolated_der, axis=1)*x_axis_epi[1]
    
    int_fluo_sums_cf_epi = np.sum(Uptake(x_axis_epi[:int(epidermis_size_mean)], fit_A_epi, fit_B_epi, fit_C_epi))*x_axis_epi[1]
    int_fluo_sums_cf_epi_err_A = int_fluo_sums_cf_epi - np.sum(Uptake(x_axis_epi[:int(epidermis_size_mean)], fit_A_epi+fit_A_err_epi, fit_B_epi, fit_C_epi))*x_axis_epi[1]
    int_fluo_sums_cf_epi_err_B = int_fluo_sums_cf_epi - np.sum(Uptake(x_axis_epi[:int(epidermis_size_mean)], fit_A_epi, fit_B_epi+fit_B_err_epi, fit_C_epi))*x_axis_epi[1]
    int_fluo_sums_cf_epi_err_C = int_fluo_sums_cf_epi - np.sum(Uptake(x_axis_epi[:int(epidermis_size_mean)], fit_A_epi, fit_B_epi, fit_C_epi+fit_C_err_epi))*x_axis_epi[1]
    int_fluo_sums_cf_epi_err = np.sqrt(int_fluo_sums_cf_epi_err_A**2 + int_fluo_sums_cf_epi_err_B**2 + int_fluo_sums_cf_epi_err_C**2)                             
     
    int_fluo_sums_cf_der = np.sum(Uptake(x_axis_der-x_axis_der[0], fit_A_der, fit_B_der, fit_C_der))*x_axis_epi[1]
    int_fluo_sums_cf_der_err_A = int_fluo_sums_cf_der - np.sum(Uptake(x_axis_der-x_axis_der[0], fit_A_der+fit_A_err_der, fit_B_der, fit_C_der))*x_axis_epi[1]
    int_fluo_sums_cf_der_err_B = int_fluo_sums_cf_der - np.sum(Uptake(x_axis_der-x_axis_der[0], fit_A_der, fit_B_der+fit_B_err_der, fit_C_der))*x_axis_epi[1]
    int_fluo_sums_cf_der_err_C = int_fluo_sums_cf_der - np.sum(Uptake(x_axis_der-x_axis_der[0], fit_A_der, fit_B_der, fit_C_der+fit_C_err_der))*x_axis_epi[1]
    int_fluo_sums_cf_der_err = np.sqrt(int_fluo_sums_cf_der_err_A**2 + int_fluo_sums_cf_der_err_B**2 + int_fluo_sums_cf_der_err_C**2)   
    
    integrated_fluorescence_epi.append(np.mean(int_fluo_sums_epi))
    integrated_fluorescence_err_epi.append(np.std(int_fluo_sums_epi))
    
    integrated_fluorescence_der.append(np.mean(int_fluo_sums_der))
    integrated_fluorescence_err_der.append(np.std(int_fluo_sums_der))
    
    integrated_fluorescence_cf_epi.append(int_fluo_sums_cf_epi)
    integrated_fluorescence_cf_err_epi.append(int_fluo_sums_cf_epi_err)
    
    integrated_fluorescence_cf_der.append(int_fluo_sums_cf_der)
    integrated_fluorescence_cf_err_der.append(int_fluo_sums_cf_der_err)
    

################################## Plot average for time point ###############             
    plt.plot(x_axis_epi, Mean_epi, label='Oral Epithelial')
    plt.fill_between(x=x_axis_epi, y1=Mean_epi-Std_epi, y2=Mean_epi+Std_epi, alpha=0.2)
    plt.plot(x_axis_der, Mean_der, label='Lamina Propria')
    plt.fill_between(x=x_axis_der, y1=Mean_der-Std_der, y2=Mean_der+Std_der, alpha=0.2)
    plt.plot(x_axis_epi[10:], Uptake( x_axis_epi[10:], fit_A_epi, fit_B_epi, fit_C_epi), color='black', linestyle='--')
    plt.plot(x_axis_der, Uptake( x_axis_der, fit_A_der, fit_B_der, fit_C_der), color='black', linestyle='--')
    plt.xlim(-10,2000)
    plt.xlabel('Depth / um')
    plt.ylabel('RFU')
    plt.title(folder+'LightOx129 Depth Penetration in Pig Tongue')
    plt.savefig(output_path+'//'+folder+' LightOx129 Depth Penetration in Pig Tongue.png')
    plt.show()

################################## Plot Repeats for time point ###############    
    for i in range(j+1):
        plt.plot(df.iloc[:,6*i], df.iloc[:,6*i+1], label='Repeat_'+str(i+1))
        plt.plot(df.iloc[:,6*i+3], df.iloc[:,6*i+4], label='Repeat_'+str(i+1))
    plt.xlim(-10,2000)
    plt.xlabel('Depth / um')
    plt.ylabel('RFU')
    plt.title(folder+'LightOx129 Depth Penetration in Pig Tongue')
    plt.legend()
    plt.savefig(output_path+'//'+folder+' Repeats LightOx129 Depth Penetration in Pig Tongue.png')
    plt.show()
    
################################## Save results ###############################
Excelwriter = pd.ExcelWriter(output_path+'\\Results.xlsx',engine="xlsxwriter")
for i, df in enumerate (df_list):    
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
    
Excelwriter.save()

################################## Plot average for each time point ###########
for i, df in enumerate(df_list):
    plt.plot(df['Depth Epi'], df['Mean Epi'], label='Oral Epithelial '+sheet_names[i])
    plt.fill_between(x=df['Depth Epi'], y1=df['Mean Epi']-df['Std Epi'], y2=df['Mean Epi']+df['Std Epi'], alpha=0.2)
    plt.plot(df['Depth Der'], df['Mean Der'], label='Lamina Propria '+sheet_names[i])
    plt.fill_between(x=df['Depth Der'], y1=df['Mean Der']-df['Std Der'], y2=df['Mean Der']+df['Std Der'], alpha=0.2)
plt.xlim(-10,2000)
plt.xlabel('Depth / um')
plt.ylabel('RFU')
plt.legend()
plt.title('LightOx129 Depth Penetration in Pig Tongue')
plt.savefig(output_path+'//LightOx129 Depth Penetration in Pig Tongue.png')
plt.show()

################################## Plot average for each time point w/ model curve###########
for i, df in enumerate(df_list):
    plt.plot(df['Depth Epi'], df['Mean Epi'], label='Oral Epithelial '+sheet_names[i])
    plt.fill_between(x=df['Depth Epi'], y1=df['Mean Epi']-df['Std Epi'], y2=df['Mean Epi']+df['Std Epi'], alpha=0.2)
    depth = df['Depth Epi'][10:]
    model = Uptake(depth, df['fit_A_Epi'][0], df['fit_B_Epi'][0], df['fit_C_Epi'][0])
    
    plt.plot(depth, model, color='black', linestyle='--')
    plt.plot(df['Depth Der'], df['Mean Der'], label='Lamina Propria '+sheet_names[i])
    plt.fill_between(x=df['Depth Der'], y1=df['Mean Der']-df['Std Der'], y2=df['Mean Der']+df['Std Der'], alpha=0.2)
    depth = df['Depth Der']
    model = Uptake(depth-depth[0], df['fit_A_Der'][0], df['fit_B_Der'][0], df['fit_C_Der'][0])
    plt.plot(depth, model, color='black', linestyle='--')
plt.xlim(-10,2000)
plt.xlabel('Depth / um')
plt.ylabel('RFU')
plt.legend()
plt.title('LightOx129 Depth Penetration in Pig Tongue')
plt.savefig(output_path+'//with models LightOx129 Depth Penetration in Pig Tongue.png')
plt.show()

################################## Plot average for each time point log #######
for i, df in enumerate(df_list):
    plt.plot(df['Depth Epi'], df['Mean Epi'], label='Oral Epithelial '+sheet_names[i])
    plt.fill_between(x=df['Depth Epi'], y1=df['Mean Epi']-df['Std Epi'], y2=df['Mean Epi']+df['Std Epi'], alpha=0.2)
    plt.plot(df['Depth Der'], df['Mean Der'], label='Lamina Propria '+sheet_names[i])
    plt.fill_between(x=df['Depth Der'], y1=df['Mean Der']-df['Std Der'], y2=df['Mean Der']+df['Std Der'], alpha=0.2)
plt.xlim(-10,2000)
plt.ylim(10,100000)
plt.xlabel('Depth / um')
plt.ylabel('RFU')
plt.yscale('log')
plt.legend()
plt.title('LightOx129 Depth Penetration in Pig Tongue')
plt.savefig(output_path+'//LOG LightOx129 Depth Penetration in Pig Tongue.png')
plt.show()

####################### Plot average for each time point log w/ model curve ###
for i, df in enumerate(df_list):
    plt.plot(df['Depth Epi'], df['Mean Epi'], label='Oral Epithelial '+sheet_names[i])
    plt.fill_between(x=df['Depth Epi'], y1=df['Mean Epi']-df['Std Epi'], y2=df['Mean Epi']+df['Std Epi'], alpha=0.2)
    depth = df['Depth Epi'][10:]
    model = Uptake(depth, df['fit_A_Epi'][0], df['fit_B_Epi'][0], df['fit_C_Epi'][0])
    
    plt.plot(depth, model, color='black', linestyle='--')
    plt.plot(df['Depth Der'], df['Mean Der'], label='Lamina Propria '+sheet_names[i])
    plt.fill_between(x=df['Depth Der'], y1=df['Mean Der']-df['Std Der'], y2=df['Mean Der']+df['Std Der'], alpha=0.2)
    depth = df['Depth Der']
    model = Uptake(depth-depth[0], df['fit_A_Der'][0], df['fit_B_Der'][0], df['fit_C_Der'][0])
    plt.plot(depth, model, color='black', linestyle='--')
plt.xlim(-10,2000)
plt.ylim(10,100000)
plt.xlabel('Depth / um')
plt.ylabel('RFU')
plt.yscale('log')
plt.legend()
plt.title('LightOx129 Depth Penetration in Pig Tongue')
plt.savefig(output_path+'//LOG with models LightOx129 Depth Penetration in Pig Tongue.png')
plt.show()

####################### Plot integrated fluorescence each time point ###########
x_labels = np.arange(len(sheet_names))
width = 0.35
fig, ax = plt.subplots()

rects1 = ax.bar(x_labels - width/2, height = integrated_fluorescence_epi, yerr = integrated_fluorescence_err_epi, label='Oral Epithelial', edgecolor='black', color='blue', width=0.35, alpha=0.5, ecolor='black', capsize=10)
rects2 = ax.bar(x_labels + width/2, height = integrated_fluorescence_der, yerr = integrated_fluorescence_err_der, label='Lamina Propria', edgecolor='black', color='orange', width=0.35, alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Integrated Fluorescence')
ax.set_title('LightOx129 Depth Penetration in Pig Tongue')
ax.set_xticks(x_labels, sheet_names)
ax.legend()
plt.savefig(output_path+'//Integrated Fluorescence.png')
plt.show()

####################### Plot integrated fluorescence each time point ###########
x_labels = np.arange(len(sheet_names))
width = 0.35
integrated_fluorescence_epi = []
integrated_fluorescence_err_epi = []
integrated_fluorescence_der = []
integrated_fluorescence_err_der = []
fig, ax = plt.subplots()
for i, df in enumerate(df_list):
    integrated_fluorescence_epi.append(np.sum(df['Mean Epi'])*df['Depth Epi'][1])
    integrated_fluorescence_err_epi.append(np.sum(df['Std Epi'])*df['Depth Epi'][1])
    integrated_fluorescence_der.append(np.sum(df['Mean Der'])*(df['Depth Der'][1] - df['Depth Der'][0]))
    integrated_fluorescence_err_der.append(np.sum(df['Std Der'])*(df['Depth Der'][1] - df['Depth Der'][0]))
rects1 = ax.bar(x_labels - width/2, height = integrated_fluorescence_epi, yerr = integrated_fluorescence_err_epi, label='Oral Epithelial', edgecolor='black', color='blue', width=0.35, alpha=0.5, ecolor='black', capsize=10)
rects2 = ax.bar(x_labels + width/2, height = integrated_fluorescence_der, yerr = integrated_fluorescence_err_der, label='Lamina Propria', edgecolor='black', color='orange', width=0.35, alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Integrated Fluorescence')
ax.set_title('LightOx129 Depth Penetration in Pig Tongue')
ax.set_xticks(x_labels, sheet_names)
ax.legend()
plt.savefig(output_path+'//Integrated Fluorescence from averages.png')
plt.show()

####################### Plot integrated fluorescence from the model each time point ###########
x_labels = np.arange(len(sheet_names))
width = 0.35
fig, ax = plt.subplots()
rects1 = ax.bar(x_labels - width/2, height = integrated_fluorescence_cf_epi, yerr = integrated_fluorescence_cf_err_epi, label='Oral Epithelial', edgecolor='black', color='blue', width=0.35, alpha=0.5, ecolor='black', capsize=10)
rects2 = ax.bar(x_labels + width/2, height = integrated_fluorescence_cf_der, yerr = integrated_fluorescence_cf_err_der, label='Lamina Propria', edgecolor='black', color='orange', width=0.35, alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Integrated Fluorescence')
ax.set_title('LightOx129 Depth Penetration in Pig Tongue (model)')
ax.set_xticks(x_labels, sheet_names)
ax.legend()
plt.savefig(output_path+'//Integrated Fluorescence (model).png')
plt.show()

################################## Save results ###############################
Tau_df = pd.DataFrame(index=sheet_names, columns={'tau epi', 'tau_err epi', 'tau der', 'tau_err der'})
taus_epi = []
taus_err_epi = []
taus_der = []
taus_err_der = []
for i, df in enumerate (df_list):    
    tau_epi = 1 / df['fit_B_Epi'][0]
    tau_err_epi = (1 / (df['fit_B_Epi'][0] - df['fit_B_err_Epi'][0])) - tau_epi
    
    tau_der = 1 / df['fit_B_Der'][0]
    tau_err_der = (1 / (df['fit_B_Der'][0] - df['fit_B_err_Der'][0])) - tau_der
    
    taus_epi.append(tau_epi)
    taus_err_epi.append(tau_err_epi)
    
    taus_der.append(tau_der)
    taus_err_der.append(tau_err_der)
    
Tau_df['tau epi'] = taus_epi
Tau_df['tau_err epi'] = taus_err_epi
Tau_df['tau der'] = taus_der
Tau_df['tau_err der'] = taus_err_der
Tau_df.to_excel(output_path+'\\Depths.xlsx', index=False)
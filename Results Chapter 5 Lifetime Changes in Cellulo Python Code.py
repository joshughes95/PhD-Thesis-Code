# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 11:40:55 2021

@author: tscw79
"""

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from tkinter.filedialog import askdirectory
import pandas as pd
import os

#Calculate Wavelengths

m=0.6471
c=383.44

x_axis=np.linspace(0,511,512)
x_axis_correct=c+m*x_axis

#first bit imports .mat data
#it also compresses cubes into 2D matrices for processing
#done both for Intensity and Lifetime Data

input_path = askdirectory(title='Select Input Directory')
output_directory = askdirectory(title='Select Output Directory')

mat_fname=input_path+'/LightOx129_5/Lifetime_Data/LifetimeImageData.mat'
mat_contents = sio.loadmat(mat_fname,squeeze_me=True)
lifetimes=mat_contents['lifetimeImageData']

mat_fname=input_path+'/LightOx129_5/Lifetime_Data/LifetimeAlphaData.mat'
mat_contents2 = sio.loadmat(mat_fname,squeeze_me=True)
intensities=mat_contents2['lifetimeAlphaData']

#Arrays
ave_lif_ROI = []
std_lif_ROI = []
int_sum_ROI = []

# Thresholding algorithm

def thres_finder(img, thres, delta_T=1.0):
    
    # Step-2: Divide the images in two parts
    x_low, y_low = np.where(img<=thres)
    x_high, y_high = np.where(img>thres)
    
    # Step-3: Find the mean of two parts
    mean_low = np.mean(img[x_low,y_low])
    mean_high = np.mean(img[x_high,y_high])
    
    # Step-4: Calculate the new threshold
    new_thres = (mean_low + mean_high)/2
    
    # Step-5: Stopping criteria, otherwise iterate
    if abs(new_thres-thres)< delta_T:
        return new_thres
    else:
        return thres_finder(img, thres=new_thres,delta_T=10)
    
# mat_fname=input_path+'/SSC4_Lo133_3/Lifetime_Data/LifetimeImageData.mat'
# mat_contents = sio.loadmat(mat_fname,squeeze_me=True)
# lifetimes=mat_contents['lifetimeImageData']

# mat_fname=input_path+'/SSC4_Lo133_3/Lifetime_Data/LifetimeAlphaData.mat'
# mat_contents2 = sio.loadmat(mat_fname,squeeze_me=True)
# intensities=mat_contents2['lifetimeAlphaData']

##################### Proocess eachrepeat in turn #############################

for h, folder in enumerate(os.listdir(input_path)):
    
    mat_fname=input_path+'/'+folder+'/Lifetime_Data/LifetimeImageData.mat'
    mat_contents = sio.loadmat(mat_fname,squeeze_me=True)
    lifetimes=mat_contents['lifetimeImageData']

    mat_fname=input_path+'/'+folder+'/Lifetime_Data/LifetimeAlphaData.mat'
    mat_contents2 = sio.loadmat(mat_fname,squeeze_me=True)
    intensities=mat_contents2['lifetimeAlphaData']

    #Produce Intensity Image
    
    sum_int=np.sum(intensities,axis=0)
    plt.figure()
    plt.imshow(sum_int)
    plt.colorbar(label='Counts')
    plt.title('Intensity')
    plt.axis('off')
    #plt.savefig(output_directory+'/'+folder+'/Summed Intensity.png', bbox_inches='tight')
    plt.show()
    
    #Intensity V Wavelength
    
    intensity_v_wave=np.sum(intensities,axis=(1,2))
    plt.figure()
    plt.plot(x_axis_correct, intensity_v_wave)
    plt.xlabel('Wavelength / nm')
    plt.ylabel('Intensity')
    #plt.savefig(output_directory+'/'+folder+'/Intensity v Wavelength.png', bbox_inches='tight')
    plt.show()
    
    # Threshold

    thres = np.average(sum_int)
    
    threshold = thres_finder(sum_int, thres)
    
    print(threshold)
    
    #Average Lifetime
    
    ave_lif=np.average(lifetimes,axis=0)
    plt.imshow(ave_lif)
    plt.colorbar(label='Lifetime (ns)')
    plt.clim(1,4)
    plt.title('Average Lifetime')
    #plt.savefig(output_directory+'/Average Lifetime.png', bbox_inches='tight')
    plt.show()
    
    #Average Lifetime, filters out lifetime = 0 spectral channels
    
    ave_lif_no_correction=np.average(lifetimes,axis=0)
    
    lifetime_image=np.zeros((len(lifetimes[0]),len(lifetimes[0][0])))
    for i in range(len(lifetimes[0])):
        for k in range(len(lifetimes[0][0])):
            lifetime_corrected = lifetimes[200:201,i,k]
            #lifetime_corrected = lifetime_corrected[lifetime_corrected != 0]
            lifetime_corrected = lifetime_corrected[lifetime_corrected > 0]
            lifetime_image[i][k]=np.average(lifetime_corrected)
    
    plt.figure()
    norm_intensity = sum_int/np.max(sum_int[:,4:])
    norm_intensity[norm_intensity < 0] = 0
    norm_intensity[norm_intensity > 1] = 1
    im1 = plt.imshow(np.ones((256,256)), cmap='Greys')
    plt.clim(0,1)
    im2 = plt.imshow(np.round(lifetime_image, 1), alpha=norm_intensity, cmap='viridis')
    plt.colorbar(label='Lifetime (ns)', shrink=0.5, aspect=7.5)
    plt.title('Average Lifetime')
    plt.clim(1.0,2.5)
    plt.axis('off') 
    #plt.savefig(output_directory+'/'+folder+'/Average Lifetime.png', bbox_inches='tight')
    plt.show()
    
    #FLIM FIGURE
    
    fig = plt.figure() # This creates the figure. 
    im1 = plt.imshow(norm_intensity, cmap='gray') # This shows the first image as a grayscale image.
    #plt.axis('off') # This removes the axes around the image.
    #plt.clim(200,800) # This clips the values within the first image to the 200-800 range. It affects the values of the colorbar.
    #plt.colorbar(shrink=0.3, aspect=5.0) # This draws a colorbar that 0.3x the size of the image and is a little bit thicker (aspect=5.0) than the default one.
    im2 = plt.imshow(lifetime_image, cmap='viridis', alpha=0.5) # This shows the second image using a 'jet' colorspace.
    plt.axis('off') # Removes the image axes.
    plt.clim(1.5,2.0) # Clips the image values to the range 2.0-3.0.
    plt.colorbar(label='Lifetime (ns)', shrink=0.5, aspect=7.5) # Draws the colorbar for the second image.
    plt.title('FLIM')
    plt.show() # Shows the figure.
    
    #Calculate and Plot lifetime vs wavelength above threshold. Run this section twice
    
    ave_lif=np.zeros((len(lifetimes)))
    std_lif=np.zeros((len(lifetimes)))
    thresholded_image = np.zeros((256, 256))
    for i in range(0,len(lifetimes)):
        plane_lif=lifetimes[i]
        above_threshold = []
        for j in range(len(plane_lif)):
            for k in range(len(plane_lif[0])):
                if sum_int[j][k] >= threshold:
                    thresholded_image[j][k]=1
                    if plane_lif[j][k] > 0:
                        above_threshold.append(plane_lif[j][k])
                
                else:
                    thresholded_image[j][k]=0
        ave_lif[i]=np.mean(above_threshold)
        std_lif[i]=np.std(above_threshold)
    
    ave_lif=np.asarray(ave_lif)
    std_lif=np.asarray(std_lif)  
    
    #ave_lif=np.average(lifetimes,axis=(1,2))
    #std_lif=np.std(lifetimes,axis=(1,2))
    
    #Plot Mask

    plt.imshow(thresholded_image, cmap='Greys')
    #plt.colorbar(label='Counts')
    plt.title('Mask')
    #plt.savefig(output_directory+'/Thresholded Image.png', bbox_inches='tight')
    plt.show()
    
    # Plot lifetime and Intensity
    
    plt.figure()
    fig,ax = plt.subplots()
    ax.plot(x_axis_correct,ave_lif, label='lifetime')
    ax.fill_between(x=x_axis_correct, y1=ave_lif-std_lif, y2=ave_lif+std_lif, alpha=0.2)
    #ax.set_ylim(0,4)
    ax.set_xlim(400, 700)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Lifetime (ns)')
    
    ax2=ax.twinx()
    ax2.plot(x_axis_correct, np.sum(intensities,axis=(1,2)), color='red', label='Intensity')
    ax2.set_ylabel('Summed Counts')
    
    fig.legend()
    #fig.savefig(output_directory+'/'+folder+'/Lifetime from Thresholded Image.png', bbox_inches='tight')
    fig.show()
    
    #Enter data into arrays
    ave_lif_ROI.append(ave_lif)
    std_lif_ROI.append(std_lif)
    int_sum_ROI.append(np.sum(intensities,axis=(1,2)))

########################## Lifetime ###########################################

ave_lif_all = np.nanmean(ave_lif_ROI, axis=0)

std_lif_all = [None]*512
numbers_to_std = np.transpose(std_lif_ROI)
numbers_to_std = numbers_to_std**2
for i in range(512):
    std_lif_all[i] = np.sqrt(np.nansum(numbers_to_std[i]))/np.count_nonzero(~np.isnan(numbers_to_std[i]))
    
ave_lif_all = np.array(ave_lif_all)
std_lif_all = np.array(std_lif_all)

###############################  Intensity ####################################

values_to_average = []
for i in range(len(int_sum_ROI)):
    values_to_average.append(int_sum_ROI[i]/np.max(int_sum_ROI[i]))
    int_sum_all = np.nanmean(values_to_average, axis=0)
    #int_sum_all = np.average((int_sum_ROI[0]/np.max(int_sum_ROI[0]), int_sum_ROI[1]/np.max(int_sum_ROI[1]), int_sum_ROI[2/np.max(int_sum_ROI[2]), int_sum_ROI[3]/np.max(int_sum_ROI[3])), axis=0)
    std_int_all = np.nanstd(values_to_average, axis=0)
    #std_int_all = np.std((int_sum_ROI[0]/np.max(int_sum_ROI[0]), int_sum_ROI[1]/np.max(int_sum_ROI[1]), int_sum_ROI[2]/np.max(int_sum_ROI[2]), int_sum_ROI[3]/np.max(int_sum_ROI[3])), axis=0)



############################### Graphs ########################################

#Repeats overlayed
colours = ['blue', 'orange', 'green', 'red', 'purple']
fig,ax = plt.subplots()
for i in range(len(ave_lif_ROI)):
    ax.plot(x_axis_correct, ave_lif_ROI[i], color=colours[i], label='ROI '+str(i+1))

ax.set_ylim(0,4)
ax.set_xlim(400, 700)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Lifetime (ns)')

ax2=ax.twinx()
for i in range(len(ave_lif_ROI)):
    ax2.plot(x_axis_correct, values_to_average[i], color=colours[i])

ax2.set_ylabel('Normalised Intensity')

fig.legend()
fig.savefig(output_directory+'/Lifetimes Overlayed.png', bbox_inches='tight')
fig.show()

#averaged plots
fig,ax = plt.subplots()
ax.plot(x_axis_correct,ave_lif_all, label='lifetime')
ax.fill_between(x=x_axis_correct, y1=ave_lif_all-std_lif_all, y2=ave_lif_all+std_lif_all, alpha=0.2)
ax.set_ylim(0,4)
ax.set_xlim(400, 700)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Lifetime (ns)')

ax2=ax.twinx()
ax2.plot(x_axis_correct, int_sum_all, color='red', label='Intensity')
ax2.fill_between(x=x_axis_correct, y1=int_sum_all-std_int_all, y2=int_sum_all+std_int_all, color='red', alpha=0.2)
ax2.set_ylabel('Normalised Intensity')

fig.legend()
fig.savefig(output_directory+'/Averaged Lifetime.png', bbox_inches='tight')
fig.show()

################################### Data Export ###############################

# Export to Excel
        
Excelwriter = pd.ExcelWriter(output_directory+'\\Results.xlsx',engine="xlsxwriter")

cols = 'Wavelength', 'Lifetime', 'Lifetime Error', 'Intensity', 'Intensity Error', 'ROI 1', 'ROI 2', 'ROI 3', 'ROI 4', 'ROI 5'

df = pd.DataFrame(data={cols[0]: x_axis_correct, cols[1]: ave_lif_all, cols[2]: std_lif_all, cols[3]: int_sum_all, cols[4]: std_int_all, cols[5]: ave_lif_ROI[0], cols[6]: ave_lif_ROI[1], cols[7]: ave_lif_ROI[2], cols[8]: ave_lif_ROI[3], cols[9]: ave_lif_ROI[4]})

df.to_excel(Excelwriter,  index=False)
Excelwriter.save()
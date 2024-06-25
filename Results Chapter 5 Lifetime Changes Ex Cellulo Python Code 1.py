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

#first bit imports .mat data
#it also compresses cubes into 2D matrices for processing
#done both for Intensity and Lifetime Data

input_path = askdirectory(title='Select Input Directory')
output_directory = askdirectory(title='Select Output Directory')

mat_fname=input_path+'/LifetimeImageData.mat'
mat_contents = sio.loadmat(mat_fname,squeeze_me=True)
lifetimes=mat_contents['lifetimeImageData']


mat_fname=input_path+'/LifetimeAlphaData.mat'
mat_contents2 = sio.loadmat(mat_fname,squeeze_me=True)
intensities=mat_contents2['lifetimeAlphaData']

#Calculate Wavelengths

m=0.5195
c=413.65

x_axis=np.linspace(0,511,512)
x_axis_correct=c+m*x_axis

#Sum Intensities

sum_int=np.sum(intensities,axis=0)
plt.imshow(sum_int)
plt.colorbar(label='Counts')
plt.title('Intensity')
plt.savefig(output_directory+'/Summed Intensity.png', bbox_inches='tight')
plt.show()

#Average Lifetime, filters out lifetime = 0 spectral channels

ave_lif_no_correction=np.average(lifetimes,axis=0)

lifetime_image=np.zeros((len(lifetimes[0]),len(lifetimes[0][0])))
for i in range(len(lifetimes[0])):
    for k in range(len(lifetimes[0][0])):
        lifetime_corrected = lifetimes[200:201,i,k]
        lifetime_corrected = lifetime_corrected[lifetime_corrected != 0]
        lifetime_image[i][k]=np.average(lifetime_corrected)

plt.imshow(np.round(lifetime_image, 1))
plt.colorbar(label='Lifetime (ns)')
plt.title('Average Lifetime')
plt.savefig(output_directory+'/Average Lifetime.png', bbox_inches='tight')
plt.show()  

#Calculate and Plot lifetime vs wavelength above threshold. Run this section twice

ave_lif=np.zeros((len(lifetimes)))
std_lif=np.zeros((len(lifetimes)))
for i in range(0,len(lifetimes)):
    plane_lif=lifetimes[i]
    above_threshold = []
    for j in range(len(plane_lif)):
        for k in range(len(plane_lif[0])):
            if plane_lif[j][k] > 0:
                above_threshold.append(plane_lif[j][k])
                
    ave_lif[i]=np.mean(above_threshold)
    std_lif[i]=np.std(above_threshold)
    
ave_lif=np.asarray(ave_lif)
std_lif=np.asarray(std_lif)   

fig,ax = plt.subplots()
ax.plot(x_axis_correct,ave_lif, label='lifetime')
ax.fill_between(x=x_axis_correct, y1=ave_lif-std_lif, y2=ave_lif+std_lif, alpha=0.2)
ax.set_ylim(0.5,2)
ax.set_xlim(400, 700)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Lifetime (ns)')

ax2=ax.twinx()
ax2.plot(x_axis_correct, np.sum(intensities,axis=(1,2)), color='red', label='Intensity')
ax2.set_ylabel('Summed Counts')

fig.legend()
fig.savefig(output_directory+'/Lifetime from Thresholded Image.png', bbox_inches='tight')
fig.show()

#Average across repeat. Run code up to here, then for each repeat run the next relevant section. Once repeats have been done, run the rest.

#Repeat 1
ave_lif_ROI1 = ave_lif
std_lif_ROI1 = std_lif
int_sum_ROI1 = np.sum(intensities,axis=(1,2))

#Repeat 2
ave_lif_ROI2 = ave_lif
std_lif_ROI2 = std_lif
int_sum_ROI2 = np.sum(intensities,axis=(1,2))

#Repeat 3
ave_lif_ROI3 = ave_lif
std_lif_ROI3 = std_lif
int_sum_ROI3 = np.sum(intensities,axis=(1,2))

#Repeat 4
ave_lif_ROI4 = ave_lif
std_lif_ROI4 = std_lif
int_sum_ROI4 = np.sum(intensities,axis=(1,2))

ave_lif_all = [None] * 512
std_lif_all = [None] * 512
for i in range(len(x_axis_correct)):
    numbers_to_average = np.array((ave_lif_ROI1[i], ave_lif_ROI2[i], ave_lif_ROI3[i], ave_lif_ROI4[i]))
    #numbers_to_average = np.array((ave_lif_ROI2[i], ave_lif_ROI3[i], ave_lif_ROI4[i]))
    #numbers_to_average = np.array((ave_lif_ROI1[i],  ave_lif_ROI3[i], ave_lif_ROI4[i]))
    numbers_to_average = numbers_to_average[numbers_to_average>0]
    ave_lif_all[i] = np.average(numbers_to_average)


    numbers_to_std = np.array((std_lif_ROI1[i], std_lif_ROI2[i], std_lif_ROI3[i], std_lif_ROI4[i]))
    #numbers_to_std = np.array((std_lif_ROI2[i], std_lif_ROI3[i], std_lif_ROI4[i]))
    #numbers_to_std = np.array((std_lif_ROI1[i], std_lif_ROI3[i],  std_lif_ROI4[i]))
    numbers_to_std = numbers_to_std[numbers_to_std>0]
    std_lif_all[i] = np.linalg.norm(numbers_to_std) / len(numbers_to_std)
    
ave_lif_all = np.array(ave_lif_all)
std_lif_all = np.array(std_lif_all)

int_sum_all = np.average((int_sum_ROI1/np.max(int_sum_ROI1), int_sum_ROI2/np.max(int_sum_ROI2), int_sum_ROI3/np.max(int_sum_ROI3), int_sum_ROI4/np.max(int_sum_ROI4)), axis=0)
std_int_all = np.std((int_sum_ROI1/np.max(int_sum_ROI1), int_sum_ROI2/np.max(int_sum_ROI2), int_sum_ROI3/np.max(int_sum_ROI3), int_sum_ROI4/np.max(int_sum_ROI4)), axis=0)

#ave_lif_all = np.average((ave_lif_ROI1, ave_lif_ROI4), axis=0)
#std_lif_all = np.sqrt(std_lif_ROI1**2 + std_lif_ROI4**2) / 2
#int_sum_all = np.average((int_sum_ROI1/np.max(int_sum_ROI1), int_sum_ROI4/np.max(int_sum_ROI4)), axis=0)
#std_int_all = np.std((int_sum_ROI1/np.max(int_sum_ROI1), int_sum_ROI4/np.max(int_sum_ROI4)), axis=0)



#Repeats overlayed
fig,ax = plt.subplots()
ax.plot(x_axis_correct,ave_lif_ROI1, color = 'blue', label='ROI 1')
ax.plot(x_axis_correct,ave_lif_ROI2, color = 'orange', label='ROI 2')
ax.plot(x_axis_correct,ave_lif_ROI3, color = 'green', label='ROI 3')
ax.plot(x_axis_correct,ave_lif_ROI4, color = 'red', label='ROI 4')

ax.set_ylim(0.5,2)
ax.set_xlim(400, 700)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Lifetime (ns)')

ax2=ax.twinx()
ax2.plot(x_axis_correct, int_sum_ROI1/np.max(int_sum_ROI1), color = 'blue')
ax2.plot(x_axis_correct, int_sum_ROI2/np.max(int_sum_ROI2), color = 'orange')
ax2.plot(x_axis_correct, int_sum_ROI3/np.max(int_sum_ROI3), color = 'green')
ax2.plot(x_axis_correct, int_sum_ROI4/np.max(int_sum_ROI4), color = 'red')
ax2.set_ylabel('Normalised Intensity')

fig.legend()
fig.savefig(output_directory+'/Lifetimes Overlayed.png', bbox_inches='tight')
fig.show()

#averaged plots
fig,ax = plt.subplots()
ax.plot(x_axis_correct,ave_lif_all, label='lifetime')
ax.fill_between(x=x_axis_correct, y1=ave_lif_all-std_lif_all, y2=ave_lif_all+std_lif_all, alpha=0.2)
ax.set_ylim(0.5,2)
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

# Export to Excel
        
Excelwriter = pd.ExcelWriter(output_directory+'\\Results.xlsx',engine="xlsxwriter")

cols = 'Wavelength', 'Lifetime', 'Lifetime Error', 'Intensity', 'Intensity Error', 'ROI 1', 'ROI 2', 'ROI 3', 'ROI4'

df = pd.DataFrame(data={cols[0]: x_axis_correct, cols[1]: ave_lif_all, cols[2]: std_lif_all, cols[3]: int_sum_all, cols[4]: std_int_all, cols[5]: ave_lif_ROI1, cols[6]: ave_lif_ROI2, cols[7]: ave_lif_ROI3, cols[8]: ave_lif_ROI4})

df.to_excel(Excelwriter,  index=False)
Excelwriter.save()

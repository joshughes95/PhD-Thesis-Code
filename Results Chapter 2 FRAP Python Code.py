# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:08:30 2020

@author: tscw79
"""
from __future__ import division
import numpy
import csv
import os
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': False})
from scipy.stats import ks_2samp

#Variables
control_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-09-02 FRAP\Control\Processed'
CR_only_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-09-02 FRAP\CR ONLY\Processed_405\data_processed.csv'
input_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-09-02 FRAP\BOTH\Processed'
output_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-09-02 FRAP\BOTH\Figures'
pre_bleach_frames = 5
Laser_powers = numpy.array([0,1,10,50])
number_of_repeats = 3
number_of_CR_repeats = 3

#Extract Control Data
control_data = []
for filename in os.listdir(control_directory):
    data_control = csv.reader(open(control_directory+'\\'+filename, 'r'), delimiter=",", quotechar='|')
    data_extracter = []
    for row in data_control:
        data_extracter.append(row)
    row_headings = data_extracter[0]
    data_array_control = numpy.zeros((len(data_extracter)-1,len(data_extracter[0])))
    for k in range(len(data_extracter)-1):
        data_array_control[k] = data_extracter[k+1]
    control_data.append(numpy.transpose(data_array_control)[1])

control_data_array = numpy.array(control_data)
control_data_averaged = numpy.zeros((len(control_data_array[0])))
for k in range(len(control_data_array[0])):
    control_data_averaged[k] = numpy.mean(control_data_array[:,k])
    
#Extract CR Only Data
CR_Only = csv.reader(open(CR_only_directory, 'r'), delimiter=",", quotechar='|')
data_extracter = []
for row in CR_Only:
    data_extracter.append(row)
row_headings = data_extracter[0]
time = numpy.array(numpy.transpose(numpy.array(data_extracter))[0], dtype=float)
CR_Only_array = numpy.transpose(numpy.array(data_extracter))[1:]

number_of_cr_repeats = int(len(CR_Only_array)/len(numpy.array(Laser_powers)))

#CR Only Arrays
cr_only_3D = numpy.zeros((len(numpy.array(Laser_powers)), number_of_cr_repeats, len(CR_Only_array[0])))
for k in range(len(numpy.array(Laser_powers))):
    for i in range(number_of_cr_repeats):
        cr_only_3D[k][i] = CR_Only_array[number_of_cr_repeats*k+i]

cr_only_data_averaged = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
cr_only_data_std = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
cr_only_data_stderr = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
cr_only_data_min = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
cr_only_data_max = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
for k in range(len(numpy.array(Laser_powers))):     
    for i in range(len(CR_Only_array[0])):
        cr_only_data_averaged[k][i] = numpy.mean(cr_only_3D[k,:,i])
        cr_only_data_std[k][i] = numpy.std(cr_only_3D[k,:,i])
        cr_only_data_stderr[k][i] = numpy.std(cr_only_3D[k,:,i])/numpy.sqrt(number_of_CR_repeats)
        cr_only_data_min[k][i] = numpy.min(cr_only_3D[k,:,i])
        cr_only_data_max[k][i] = numpy.max(cr_only_3D[k,:,i])

#Loop over each repeat
processed_data = []
area_data = []
for filename in os.listdir(input_directory):
    #Extract ROI and Cell area sizes
    if filename.startswith('AREA'):
        data = csv.reader(open(input_directory+'\\'+filename, 'r'), delimiter=",", quotechar='|')
        #Read CSV, convert to an array and transpose
        data_extracter = []
        for row in data:
            data_extracter.append(row)
        data_array = numpy.zeros((len(data_extracter)-1,len(data_extracter[0])))
        for i in range(len(data_extracter)-1):
            data_array[i] = data_extracter[i+1]
        area_data.append(numpy.transpose(data_array)[0])
        
for filename in os.listdir(input_directory):
    #Extract FRAP raw data from CSV
    if filename.startswith('FRAP'):
        data = csv.reader(open(input_directory+'\\'+filename, 'r'), delimiter=",", quotechar='|')
        #Read CSV, convert to an array and transpose
        data_extracter = []
        for row in data:
            data_extracter.append(row)
        data_array = numpy.zeros((len(data_extracter)-1,len(data_extracter[0])))
        for i in range(len(data_extracter)-1):
            data_array[i] = data_extracter[i+1]
        
        #Initial Arrays
        data_sorted = numpy.transpose(data_array)
        data_raw = data_sorted[1: -1]
        data_fade = data_sorted[-1]
        
        #find the correct Cell Area for the FRAP Data
        for i in range(int(len(os.listdir(input_directory))/2)+1):
            if filename.endswith(str(i+1)+'.csv'):
                ROI_area = area_data[i][0]
                cell_area = area_data[i][1]
        
        #Correct Fade and Background for bleaching ROIs
        data_fade_corrected = numpy.zeros((1,len(data_fade)))[0]
        sum_rows = [sum(x) for x in zip(*data_raw)]
        for i in range(len(data_fade)):
            data_fade_corrected[i] = (data_fade[i]*cell_area-(ROI_area*sum_rows[i]))/(cell_area-3*ROI_area)
        data_corrected = (data_raw - control_data_averaged) / (data_fade_corrected - control_data_averaged)
        
        #Normalise Data
        normalisation_constants = numpy.zeros((1,len(data_corrected)))[0]
        data_normalised = numpy.zeros((len(data_corrected),len(data_corrected[0])))
        for i in range(len(data_corrected)):
            normalisation_constants[i] += numpy.mean(data_corrected[i][1:pre_bleach_frames-1])
            data_normalised[i] = data_corrected[i]/normalisation_constants[i]
        processed_data.append(data_normalised)

data_processed = numpy.array(processed_data)
data_processed_averaged = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
data_processed_std = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
data_processed_stderr = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
data_processed_min = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
data_processed_max = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
for k in range(len(data_processed[0])):
    data_processed_averaged[k] = [numpy.mean(x) for x in zip(*data_processed[:,k])]
    data_processed_std[k] = [numpy.std(x) for x in zip(*data_processed[:,k])]
    data_processed_stderr[k] = [numpy.std(x) for x in zip(*data_processed[:,k])]/numpy.sqrt(number_of_repeats)
    data_processed_min[k] = [numpy.min(x) for x in zip(*data_processed[:,k])]
    data_processed_max[k] = [numpy.max(x) for x in zip(*data_processed[:,k])]

#Correct for CR ONLY fading
#data_processed_average = data_processed_averaged / cr_only_data_averaged
data_processed_average = data_processed_averaged

#kolmogorov-smirnov stat
D = 1.36*numpy.sqrt(2/(len(data_processed_averaged[0])-1))
print(D)
file = open(output_directory+'\\'+'Less than 75'+str(Laser_powers[i])+'.txt', 'w')
for i in range(len(Laser_powers)):
    KstestResult = ks_2samp(data_processed_averaged[i][1:], cr_only_data_averaged[i][1:])
    print(KstestResult)
    if KstestResult[0] < D:
        file.writelines('The data sets are not distinct for Laser Power {} nm. Statistic = {}, P-value = {}\n'.format(Laser_powers[i], KstestResult[0], KstestResult[1]))
    else:
        file.writelines('The data sets are distinct for Laser Power {} nm. Statistic = {}, P-value ={}\n'.format(Laser_powers[i], KstestResult[0], KstestResult[1]))

#Propagate Errors (stderr only)
data_propagated_error = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
for k in range(len(data_processed[0])):
    for i in range(len(data_processed[0][0])):
        data_propagated_error[k][i] = numpy.linalg.norm([data_processed_stderr[k][i],cr_only_data_stderr[k][i]])
    
#Add time to file
data_with_time = numpy.zeros((len(processed_data)*len(processed_data[0])+1, len(processed_data[0][0])))
data_with_time[0] = time
for k in range(len(data_processed)):
    for i in range(len(data_processed[0])):
        data_with_time[len(data_processed[0])*k+1+i] = data_processed[k][i]

#Create file with processed data
data_transposed = numpy.transpose(data_with_time)
with open(input_directory+'\\'+'data_processed.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for i in range(len(data_with_time[0])):
        writer.writerow(data_transposed[i])

# for i in range(len(data_processed_average)):
#     plt.errorbar(x=time, y=data_processed_average[i], yerr=data_propagated_error[i], color='black', ecolor='tab:blue')
#     plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
#     plt.vlines(x=(time[pre_bleach_frames]+time[pre_bleach_frames-1])/2, ymin=0.5, ymax=1.5, linestyle='--', color= 'grey')
#     plt.ylim(0.5,1.5)
#     plt.xlabel('Time/s')
#     plt.ylabel('Normalised Intensity/au')
#     plt.xticks(numpy.arange(float(min(time)), float(max(time)), 20))
#     plt.title('Laser Power = %d' %Laser_powers[i])
#     plt.savefig(output_directory+'\\'+'Laser_power_'+str(Laser_powers[i])+'.png')
#     plt.show()
    
#Figures
#std

for i in range(len(data_processed_average)):
    plt.errorbar(x=time, y=data_processed_average[i], yerr=data_processed_std[i], color='black', ecolor='tab:blue')
    #plt.plot(time, data_processed_average[i], color='black')
    #plt.fill_between(x=time, y1=data_processed_average[i] - data_processed_std[i], y2=data_processed_average[i] + data_processed_std[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=(time[pre_bleach_frames]+time[pre_bleach_frames-1])/2, ymin=0.5, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.5,1.5)
    plt.xlabel('Time/s')
    plt.ylabel('Normalised Intensity/au')
    plt.title('Laser Power = %d' %Laser_powers[i])
    plt.savefig(output_directory+'\\'+'(std) Laser_power_'+str(Laser_powers[i])+'.png')
    plt.show()
 
#std_err
for i in range(len(data_processed_average)):
    plt.errorbar(x=time, y=data_processed_average[i], yerr=data_processed_stderr[i], color='black', ecolor='tab:blue')
    #plt.plot(time, data_processed_average[i], color='black')
    #plt.fill_between(x=time, y1=data_processed_average[i] - data_processed_stderr[i], y2=data_processed_average[i] + data_processed_stderr[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=(time[pre_bleach_frames]+time[pre_bleach_frames-1])/2, ymin=0.5, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.5,1.5)
    plt.xlabel('Time/s')
    plt.ylabel('Normalised Intensity/au')
    plt.title('Laser Power = %d' %Laser_powers[i])
    plt.savefig(output_directory+'\\'+'(std_err) Laser_power_'+str(Laser_powers[i])+'.png')
    plt.show()

#ranges  
for i in range(len(data_processed_average)):
    plt.errorbar(x=time, y=data_processed_average[i], yerr=(data_processed_average[i] - data_processed_min[i], data_processed_max[i] - data_processed_average[i]), color='black', ecolor='tab:blue', label='Measured', alpha=0.6)
    plt.plot(time, data_processed_average[i], color='black')
    #plt.fill_between(x=time, y1=data_processed_min[i], y2=data_processed_max[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=(time[pre_bleach_frames]+time[pre_bleach_frames-1])/2, ymin=0.5, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.5,1.5)
    plt.xlabel('Time/s')
    plt.ylabel('Normalised Intensity/au')
    plt.title('Laser Power = %d' %Laser_powers[i])
    plt.savefig(output_directory+'\\'+'(Ranges) Laser_power_'+str(Laser_powers[i])+'.png')
    plt.show()
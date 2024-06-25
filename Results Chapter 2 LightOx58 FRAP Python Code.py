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
import scipy
from scipy.optimize import minimize
plt.rcParams.update({'font.size': 16})

#Variables
#control_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-08-27 FRAP\Control\Processed'
#CR_only_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2020-08-27 FRAP\CR ONLY\Processed_405\data_processed.csv'
input_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\LightOx58 FRAP\Processed_less_time'
output_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\LightOx58 FRAP\Figures_less_time'
pre_bleach_frames = 5
Laser_powers = numpy.array([405, 405633, 633, 0])
number_of_repeats = 3
#number_of_CR_repeats = 3

#Extract Control Data
# control_data = []
# for filename in os.listdir(control_directory):
#     data_control = csv.reader(open(control_directory+'\\'+filename, 'r'), delimiter=",", quotechar='|')
#     data_extracter = []
#     for row in data_control:
#         data_extracter.append(row)
#     row_headings = data_extracter[0]
#     data_array_control = numpy.zeros((len(data_extracter)-1,len(data_extracter[0])))
#     for k in range(len(data_extracter)-1):
#         data_array_control[k] = data_extracter[k+1]
#     control_data.append(numpy.transpose(data_array_control)[1])

# control_data_array = numpy.array(control_data)
# control_data_averaged = numpy.zeros((len(control_data_array[0])))
# for k in range(len(control_data_array[0])):
#     control_data_averaged[k] = numpy.mean(control_data_array[:,k])
    
# #Extract CR Only Data
# CR_Only = csv.reader(open(CR_only_directory, 'r'), delimiter=",", quotechar='|')
# data_extracter = []
# for row in CR_Only:
#     data_extracter.append(row)
# row_headings = data_extracter[0]
# time = numpy.array(numpy.transpose(numpy.array(data_extracter))[0], dtype=float)
# CR_Only_array = numpy.transpose(numpy.array(data_extracter))[1:]

# number_of_cr_repeats = int(len(CR_Only_array)/len(numpy.array(Laser_powers)))

# #CR Only Arrays
# cr_only_3D = numpy.zeros((len(numpy.array(Laser_powers)), number_of_cr_repeats, len(CR_Only_array[0])))
# for k in range(len(numpy.array(Laser_powers))):
#     for i in range(number_of_cr_repeats):
#         cr_only_3D[k][i] = CR_Only_array[number_of_cr_repeats*k+i]

# cr_only_data_averaged = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
# cr_only_data_std = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
# cr_only_data_stderr = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
# cr_only_data_min = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
# cr_only_data_max = numpy.zeros((len(Laser_powers), len(CR_Only_array[0])))
# for k in range(len(numpy.array(Laser_powers))):     
#     for i in range(len(CR_Only_array[0])):
#         cr_only_data_averaged[k][i] = numpy.mean(cr_only_3D[k,:,i])
#         cr_only_data_std[k][i] = numpy.std(cr_only_3D[k,:,i])
#         cr_only_data_stderr[k][i] = numpy.std(cr_only_3D[k,:,i])/numpy.sqrt(number_of_CR_repeats)
#         cr_only_data_min[k][i] = numpy.min(cr_only_3D[k,:,i])
#         cr_only_data_max[k][i] = numpy.max(cr_only_3D[k,:,i])

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
        time = numpy.array(numpy.transpose(numpy.array(data_extracter[1:]))[0], dtype=float)
        #data_array = numpy.array(data_extracter[1:], dtype=float)
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
        data_corrected = (data_raw) / (data_fade_corrected)
        
        #Normalise Data
        normalisation_constants = numpy.zeros((1,len(data_corrected)))[0]
        data_normalised = numpy.zeros((len(data_corrected),len(data_corrected[0])))
        for i in range(len(data_corrected)):
            normalisation_constants[i] += numpy.mean(data_corrected[i][1:pre_bleach_frames-1])
            data_normalised[i] = data_corrected[i]/normalisation_constants[i]
        processed_data.append(data_normalised)

data_processed = numpy.zeros((int(len(processed_data)), int(len(processed_data[0])), len(processed_data[0][0])))
for k in range(len(processed_data)):
    for i in range(len(data_processed[0])):
        for j in range(len(data_processed[0][0])):
            data_processed[k][i][j] = processed_data[k][i][j]

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

#Propagate Errors (stderr only)
# data_propagated_error = numpy.zeros((len(data_processed[0]),len(data_processed[0][0])))
# for k in range(len(data_processed[0])):
#     for i in range(len(data_processed[0][0])):
#         data_propagated_error[k][i] = numpy.linalg.norm([data_processed_stderr[k][i],cr_only_data_stderr[k][i]])


#Fit to Model


x0_linear =  numpy.array([[1,0], [1,0], [1,0], [1,0]])
x0_monoexp =  numpy.array([[1,-0.5,0.01], [1,-0.5,0.01], [1,0,0.01], [1,0,0.01]])
x0_biexp =  numpy.array([[1,-0.5,0.01,-0.5,0.01],[1,-0.5,0.01,-0.5,0.01], [1,0,0.01,0,0.01], [1,0,0.01,0,0.01]])

model_linear = []
model_mono = []
model_bi =[]

time_for_model = time[pre_bleach_frames:] - time[pre_bleach_frames]

for i in range(len(Laser_powers)):
    def chisqfunclinear(x):
        'Linear chi squared fitting function'
        a,b = x
        model = a + b*time_for_model
        chisq = numpy.sum(((data_processed_averaged[i][pre_bleach_frames:] - model)/data_processed_stderr[i][pre_bleach_frames:])**2)
        return chisq
    
    def chisqfuncmonoexp(x):
        'Linear chi squared fitting function'
        a,b,c = x
        model = a + b*numpy.exp(-c*time_for_model)
        chisq = numpy.sum(((data_processed_averaged[i][pre_bleach_frames:] - model)/data_processed_stderr[i][pre_bleach_frames:])**2)
        return chisq
    
    def chisqfuncbiexp(x):
        'Linear chi squared fitting function'
        a,b,c,d,e = x
        model = a + b*numpy.exp(-c*time_for_model) + d*numpy.exp(-e*time_for_model)
        chisq = numpy.sum(((data_processed_averaged[i][pre_bleach_frames:] - model)/data_processed_stderr[i][pre_bleach_frames:])**2)
        return chisq

    result_linear =  minimize(chisqfunclinear, x0_linear[i])
    model_linear.append(result_linear.x[0] + result_linear.x[1]*time_for_model)
    
    #a,b=result_linear.x
    
    result_mono =  minimize(chisqfuncmonoexp, x0_monoexp[i])
    model_mono.append(result_mono.x[0] + result_mono.x[1]*numpy.exp(-result_mono.x[2]*time_for_model))
    MF_mono = (result_mono.x[0] - data_processed_averaged[i][pre_bleach_frames]) / (1 - data_processed_averaged[i][pre_bleach_frames])
    #a,b=result_mono.x
    
    result_bi =  minimize(chisqfuncbiexp, x0_biexp[i])
    model_bi.append(result_bi.x[0] + result_bi.x[1]*numpy.exp(-result_bi.x[2]*time_for_model) + result_bi.x[3]*numpy.exp(-result_bi.x[4]*time_for_model))
    MF_bi = (result_bi.x[0] - data_processed_averaged[i][pre_bleach_frames]) / (1 - data_processed_averaged[i][pre_bleach_frames])
    #a,b=result_bi.x
    
    print('The Mobile Fraction from the Mono Exp Model is {} for Laser Line {} nm'.format(MF_mono, Laser_powers[i]))
    print('The Mobile Fraction from the Bi Exp Model is {} for the Laser line {} nm'.format(MF_bi, Laser_powers[i]))
    
    file = open(output_directory+'\\'+'MF_laser_line'+str(Laser_powers[i])+'.txt', 'w')
    file.writelines('The Linear model is {} + {}t for Laser Line {} nm \n'.format(result_linear.x[0], result_linear.x[1], Laser_powers[i]))
    file.writelines('The MonoExp model is {} + {}exp(-{}t) for Laser Line {} nm \n'.format(result_mono.x[0], result_mono.x[1], result_mono.x[2], Laser_powers[i]))
    file.writelines('The BiExp model is {} + {}exp(-{}t) + {}exp(-{}t) for Laser Line {} nm \n'.format(result_bi.x[0], result_bi.x[1], result_bi.x[2], result_bi.x[3], result_bi.x[4], Laser_powers[i]))
    file.writelines('The Mobile Fraction from the Mono Exp Model is {} for Laser Line {} nm \n'.format(MF_mono, Laser_powers[i]))
    file.writelines('The Mobile Fraction from the Bi Exp Model is {} for the Laser line {} nm \n'.format(MF_bi, Laser_powers[i]))

    
    plt.figure()
    chis_x = numpy.array([0, 1, 2])
    dOF_linear = len(data_processed_averaged[0]) - pre_bleach_frames - 2
    dOF_mono = len(data_processed_averaged[0]) - pre_bleach_frames - 3
    dOF_bi = len(data_processed_averaged[0]) - pre_bleach_frames - 5
    chis_y = numpy.array([chisqfunclinear(result_linear.x)/dOF_linear, chisqfuncmonoexp(result_mono.x)/dOF_mono, chisqfuncbiexp(result_bi.x)/dOF_bi])
    plt.plot(chis_x, chis_y)
    plt.title('chi square for linear, mono exp and bi exp models')
    plt.savefig(output_directory+'\\'+'Reduced Chi Square Comparisons Laser Line'+str(Laser_powers[i])+'.png')
    plt.show()
    
    corr_matrix_linear = numpy.corrcoef(data_processed_averaged[i][pre_bleach_frames:], model_linear)
    corr_matrix_mono = numpy.corrcoef(data_processed_averaged[i][pre_bleach_frames:], model_mono)
    corr_matrix_bi = numpy.corrcoef(data_processed_averaged[i][pre_bleach_frames:], model_bi)
    corr_linear = corr_matrix_linear[0,1]
    corr_mono = corr_matrix_mono[0,1]
    corr_bi = corr_matrix_bi[0,1]
    file = open(output_directory+'\\'+'RSquared_laser_line'+str(Laser_powers[i])+'.txt', 'w')
    file.writelines('The Linear Model RSquared is {} for Laser Line {} nm \n'.format(corr_linear**2, Laser_powers[i]))
    file.writelines('The Mono Exp Model RSquared is {} for Laser Line {} nm \n'.format(corr_mono**2, Laser_powers[i]))
    file.writelines('The Bi Exp Model RSquared is {} for Laser Line {} nm \n'.format(corr_bi**2, Laser_powers[i]))
    
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
    plt.errorbar(x=time, y=data_processed_average[i], yerr=data_processed_std[i], label='Measured', color='black', ecolor='tab:blue', alpha=0.4)
    plt.plot(time[pre_bleach_frames:], model_linear[i], label = 'linear')
    plt.plot(time[pre_bleach_frames:], model_mono[i], label = 'mono')
    plt.plot(time[pre_bleach_frames:], model_bi[i], label = 'bi')
    #plt.plot(time, data_processed_average[i], color='black')
    #plt.fill_between(x=time, y1=data_processed_average[i] - data_processed_std[i], y2=data_processed_average[i] + data_processed_std[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=time[pre_bleach_frames], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.vlines(x=time[pre_bleach_frames-1], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.4,1.2)
    plt.xlabel('Time/s')
    plt.ylabel('Normalised Intensity/au')
    plt.title('Laser Line = %d nm' %Laser_powers[i])
    plt.legend(title='Models', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_directory+'\\'+'(std) Laser_Line'+str(Laser_powers[i])+'.png', bbox_inches='tight')
    plt.show()
 
#std_err
for i in range(len(data_processed_average)):
    plt.errorbar(x=time, y=data_processed_average[i], yerr=data_processed_stderr[i], color='black', ecolor='tab:blue', label='Measured', alpha=0.8)
    plt.plot(time[pre_bleach_frames:], model_linear[i], label = 'Model', color = 'orange')
    #plt.plot(time[pre_bleach_frames:], model_mono[i], label = 'mono', color='orange')
    #plt.plot(time[pre_bleach_frames:], model_bi[i], label = 'bi', color='orange')
    #plt.plot(time, data_processed_average[i], color='black')
    #plt.fill_between(x=time, y1=data_processed_average[i] - data_processed_stderr[i], y2=data_processed_average[i] + data_processed_stderr[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=time[pre_bleach_frames], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.vlines(x=time[pre_bleach_frames-1], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.4,1.2)
    plt.xlabel('Time/s')
    plt.ylabel('Normalised Intensity/au')
    plt.title('Laser Line = %d nm' %Laser_powers[i])
    #plt.legend(title='Models', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_directory+'\\'+'(std_err) Laser_Line_'+str(Laser_powers[i])+'.png', bbox_inches='tight')
    plt.show()

#ranges  
for i in range(len(data_processed_average)):
    fig1 = plt.figure()
    frame1 = fig1.add_axes((0.1,0.3,0.8,0.8))
    plt.errorbar(x=time, y=data_processed_average[i], yerr=(data_processed_average[i] - data_processed_min[i], data_processed_max[i] - data_processed_average[i]), color='black', ecolor='tab:blue', label='Measured', alpha=0.6)
    plt.plot(time, data_processed_average[i], color='black', label='Measured')
    #plt.plot(time[pre_bleach_frames:], model_linear[i], label = 'linear', color='orange')
    #plt.plot(time[pre_bleach_frames:], model_mono[i], label = 'mono')
    plt.plot(time[pre_bleach_frames:], model_bi[i], label = 'bi', color='orange')
    #plt.fill_between(x=time, y1=data_processed_min[i], y2=data_processed_max[i], color='blue', alpha=0.2)
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.vlines(x=time[pre_bleach_frames], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.vlines(x=time[pre_bleach_frames-1], ymin=0, ymax=1.5, linestyle='--', color= 'grey')
    plt.ylim(0.35,1.2)
    plt.xlim(-5,190)
    plt.xlabel('Time/s')
    plt.ylabel('Normalized Intensity/au')
    plt.title('Laser Line = %d nm' %Laser_powers[i])
    
    frame2 = fig1.add_axes((0.1,0.1,0.8,0.2))
    plt.plot(time[pre_bleach_frames:], (data_processed_average[i][pre_bleach_frames:]-model_linear[i])/data_processed_stderr[i][pre_bleach_frames:], linestyle = "", marker = "o", markersize = 2)
    plt.xlabel('Time/s')
    plt.ylabel('Residuals')
    plt.xlim(-5,190)
    plt.ylim(-2.2,2.2)
    #plt.legend(title='Models', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_directory+'\\'+'(Ranges) Laser_Line'+str(Laser_powers[i])+'.png', bbox_inches='tight')
    plt.show()
    
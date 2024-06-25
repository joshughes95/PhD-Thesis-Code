# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 23:43:21 2021

@author: tscw79
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from tkinter.filedialog import askdirectory
import pandas as pd
from sklearn.linear_model import LinearRegression

input_path = askdirectory(title='Select Input Directory')
output_directory = askdirectory(title='Select Output Directory')

# define spectral window to average lifetime over, y = 413.65 + 0.5195x
wavelength_start = np.array((475, 498, 535, 576))
wavelength_end = np.array((498, 535, 576, 583))
spectral_channel_start = np.round((wavelength_start - 413.65) / 0.5195)
spectral_channel_end = np.round((wavelength_end - 413.65) / 0.5195)

for j in range(len(wavelength_start)):
    # Empty arrays for lifetime v Conc plots
    Concentrations = []
    av_lifetime_Con = []
    av_lifetime_er_Con = []
    av_lifetime_HDAC1 = []
    av_lifetime_er_HDAC1 = []
    av_lifetime_HDAC2 = []
    av_lifetime_er_HDAC2 = []
    
    # Empty arrays for Intensity v Conc plots
    max_intensity_Con = []
    max_intensity_HDAC1 = []
    max_intensity_HDAC2 = []
    
    # Empty arrays for lifetime v wavelength plots
    Control_concs = []
    HDAC1_concs = []
    HDAC2_concs =[]
    
    Control_FS = []
    HDAC1_FS = []
    HDAC2_FS =[]
    
    Control_label = []
    HDAC1_label = []
    HDAC2_label = []
    
    for filename in os.listdir(input_path):
        df = pd.read_excel(input_path+'\\'+filename, sheet_name = None)  
    
        wavelengths = df['Sheet1']['Wavelength']
        lifetime = df['Sheet1']['Lifetime']
        lifetime_error = df['Sheet1']['Lifetime Error']
        intensity = np.array(df['Sheet1']['Intensity'])
        intensity_error = df['Sheet1']['Intensity Error']
        repeat_1 = df['Sheet1']['ROI 1']
        repeat_2 = df['Sheet1']['ROI 2']
        repeat_3 = df['Sheet1']['ROI 3']
        repeat_4 = df['Sheet1']['ROI4']
        
        if 'Control' in filename:
            if '0_125uM' in filename:
                Control_FS.append(lifetime)
                Concentrations.append(0.125)
                Control_concs.append(0.125)
                Control_label.append('Control 0.125uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_Con.append(np.average(lifetime))
                av_lifetime_er_Con.append(np.std(lifetime))
                max_intensity_Con.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_25uM' in filename:
                Control_FS.append(lifetime)
                Concentrations.append(0.25)
                Control_concs.append(0.25)
                Control_label.append('Control 0.25uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_Con.append(np.average(lifetime))
                av_lifetime_er_Con.append(np.std(lifetime))
                max_intensity_Con.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_5uM' in filename:
                Control_FS.append(lifetime)
                Concentrations.append(0.5)
                Control_concs.append(0.5)
                Control_label.append('Control 0.5uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_Con.append(np.average(lifetime))
                av_lifetime_er_Con.append(np.std(lifetime))
                max_intensity_Con.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '1uM' in filename:
                Control_FS.append(lifetime)
                Concentrations.append(1.0)
                Control_concs.append(1.0)
                Control_label.append('Control 1.0uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_Con.append(np.average(lifetime))
                av_lifetime_er_Con.append(np.std(lifetime))
                max_intensity_Con.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '2uM' in filename:
                Control_FS.append(lifetime)
                Concentrations.append(2.0)
                Control_concs.append(2.0)
                Control_label.append('Control 2.0uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_Con.append(np.average(lifetime))
                av_lifetime_er_Con.append(np.std(lifetime))
                max_intensity_Con.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
                
        if 'HDAC1' in filename:
            if '0_125uM' in filename:
                HDAC1_FS.append(lifetime)
                HDAC1_concs.append(0.125)
                HDAC1_label.append('HDAC1 0.125uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC1.append(np.average(lifetime))
                av_lifetime_er_HDAC1.append(np.std(lifetime))
                max_intensity_HDAC1.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_25uM' in filename:
                HDAC1_FS.append(lifetime)
                HDAC1_concs.append(0.25)
                HDAC1_label.append('HDAC1 0.25uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC1.append(np.average(lifetime))
                av_lifetime_er_HDAC1.append(np.std(lifetime))
                max_intensity_HDAC1.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_5uM' in filename:
                HDAC1_FS.append(lifetime)
                HDAC1_concs.append(0.5)
                HDAC1_label.append('HDAC1 0.5uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC1.append(np.average(lifetime))
                av_lifetime_er_HDAC1.append(np.std(lifetime))
                max_intensity_HDAC1.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '1uM' in filename:
                HDAC1_FS.append(lifetime)
                HDAC1_concs.append(1.0)
                HDAC1_label.append('HDAC1 1uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC1.append(np.average(lifetime))
                av_lifetime_er_HDAC1.append(np.std(lifetime))
                max_intensity_HDAC1.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '2uM' in filename:
                HDAC1_FS.append(lifetime)
                HDAC1_concs.append(2.0)
                HDAC1_label.append('HDAC1 2uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC1.append(np.average(lifetime))
                av_lifetime_er_HDAC1.append(np.std(lifetime))
                max_intensity_HDAC1.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
                
        if 'HDAC2' in filename:
            if '0_125uM' in filename:
                HDAC2_FS.append(lifetime)
                HDAC2_concs.append(0.125)
                HDAC2_label.append('HDAC2 0.125uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC2.append(np.average(lifetime))
                av_lifetime_er_HDAC2.append(np.std(lifetime))
                max_intensity_HDAC2.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_25uM' in filename:
                HDAC2_FS.append(lifetime)
                HDAC2_concs.append(0.25)
                HDAC2_label.append('HDAC2 0.25uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC2.append(np.average(lifetime))
                av_lifetime_er_HDAC2.append(np.std(lifetime))
                max_intensity_HDAC2.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '0_5uM' in filename:
                HDAC2_FS.append(lifetime)
                HDAC2_concs.append(0.5)
                HDAC2_label.append('HDAC2 0.5uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC2.append(np.average(lifetime))
                av_lifetime_er_HDAC2.append(np.std(lifetime))
                max_intensity_HDAC2.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '1uM' in filename:
                HDAC2_FS.append(lifetime)
                HDAC2_concs.append(1.0)
                HDAC2_label.append('HDAC2 1uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC2.append(np.average(lifetime))
                av_lifetime_er_HDAC2.append(np.std(lifetime))
                max_intensity_HDAC2.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
            elif '2uM' in filename:
                HDAC2_FS.append(lifetime)
                HDAC2_concs.append(2.0)
                HDAC2_label.append('HDAC2 2uM')
                lifetime = lifetime[int(spectral_channel_start[j]):int(spectral_channel_end[j])]
                lifetime = lifetime[lifetime>0]
                av_lifetime_HDAC2.append(np.average(lifetime))
                av_lifetime_er_HDAC2.append(np.std(lifetime))
                max_intensity_HDAC2.append(wavelengths[int(np.array(np.where(intensity == np.max(intensity))))])
                
    ################### Convert Arrays ############################################
    
    Control_concs = np.array(Control_concs)
    HDAC1_concs = np.array(HDAC1_concs)
    HDAC2_concs = np.array(HDAC2_concs)
    Concentrations = np.array(Concentrations)
    av_lifetime_Con = np.array(av_lifetime_Con)
    av_lifetime_er_Con = np.array(av_lifetime_er_Con)
    av_lifetime_HDAC1 = np.array(av_lifetime_HDAC1)
    av_lifetime_er_HDAC1 = np.array(av_lifetime_er_HDAC1)
    av_lifetime_HDAC2 = np.array(av_lifetime_HDAC2)
    av_lifetime_er_HDAC2 = np.array(av_lifetime_er_HDAC2)
    max_intensity_Con = np.array(max_intensity_Con)
    max_intensity_HDAC1 = np.array(max_intensity_HDAC1)
    max_intensity_HDAC2 = np.array(max_intensity_HDAC2)
    
    ######################## Data Outputs #########################################
    cols = 'Concentration / uM', 'Lifetime Con / ns', 'Lifetime HDAC1 / ns', 'Lifetime HDAC2 / ns', 'LC Error / ns', 'LH1 Error / ns', 'LH2 Error / ns', 'Peak Emission Con / nm', 'Peak Emission HDAC1 / nm', 'Peak Emission HDAC2 / nm'
    
    # Export data to Excel
    
    Excelwriter = pd.ExcelWriter(output_directory+'\\Results_'+str(wavelength_start[j])+'to'+str(wavelength_end[j])+'nm.xlsx',engine="xlsxwriter")
        
    df = pd.DataFrame(data={cols[0]: Concentrations, cols[1]: av_lifetime_Con, cols[2]: av_lifetime_er_Con, cols[3]: av_lifetime_HDAC1, cols[4]: av_lifetime_er_HDAC1, cols[5]: av_lifetime_HDAC2, cols[6]: av_lifetime_er_HDAC2, cols[7]: max_intensity_Con, cols[8]: max_intensity_HDAC1, cols[9]: max_intensity_HDAC2})
    
    df.to_excel(Excelwriter,  index=False)
    Excelwriter.save()
    
    ################################# Graphs ######################################
    #Concentration v lifetime
    plt.figure()
    model = LinearRegression().fit(Control_concs.reshape((-1, 1)), av_lifetime_Con)
    coeff = model.coef_
    r_sq = model.score(Control_concs.reshape((-1, 1)), av_lifetime_Con)
    plt.errorbar(x=Control_concs,y=av_lifetime_Con, yerr=av_lifetime_er_Con, label='Control, Gradient='+str(np.round(coeff, decimals=2))+', R^2='+str(np.round(r_sq, decimals=2)), capsize=4)
    
    model = LinearRegression().fit(Control_concs.reshape((-1, 1)), av_lifetime_HDAC1)
    coeff = model.coef_
    r_sq = model.score(Control_concs.reshape((-1, 1)), av_lifetime_HDAC1)
    plt.errorbar(x=HDAC1_concs,y=av_lifetime_HDAC1, yerr=av_lifetime_er_HDAC1, label='HDAC1, Gradient='+str(np.round(coeff, decimals=2))+', R^2='+str(np.round(r_sq, decimals=2)), capsize=4)
    
    model = LinearRegression().fit(Control_concs.reshape((-1, 1)), av_lifetime_HDAC2)
    coeff = model.coef_
    r_sq = model.score(Control_concs.reshape((-1, 1)), av_lifetime_HDAC2)
    plt.errorbar(x=HDAC2_concs,y=av_lifetime_HDAC2, yerr=av_lifetime_er_HDAC2, label='HDAC2, Gradient='+str(np.round(coeff, decimals=2))+', R^2='+str(np.round(r_sq, decimals=2)), capsize=4)
    
    plt.ylim(1.0,2.0)
    plt.xlabel('LightOx Concentration [uM]')
    plt.ylabel('Lifetime (ns)')
    plt.title(str(wavelength_start[j])+'nm to '+str(wavelength_end[j])+'nm')
    plt.legend()
    
    plt.savefig(output_directory+'/Concentration v lifetime'+str(wavelength_start[j])+'to'+str(wavelength_end[j])+'nm.png', bbox_inches='tight')
    plt.show()
    
    #Intensity v lifetime
    plt.figure()
    plt.plot(Control_concs,max_intensity_Con, label='Control')
    plt.plot(Control_concs,max_intensity_HDAC1, label='HDAC1')
    plt.plot(Control_concs,max_intensity_HDAC2, label='HDAC2')
    
    plt.ylim(510,530)
    plt.xlabel('LightOx Concentration [uM]')
    plt.ylabel('Peak Emission / nm')
    plt.legend()
    
    plt.savefig(output_directory+'/Intensity v lifetime.png', bbox_inches='tight')
    plt.show()
    
    
    #overlayed plots ALL Repeats
    plt.figure()
    for i in range (len(Control_concs)):
        plt.plot(wavelengths, Control_FS[i], label = Control_label[i])
    
    for i in range (len(HDAC1_concs)):
        plt.plot(wavelengths, HDAC1_FS[i], label = HDAC1_label[i])
        
    for i in range (len(HDAC2_concs)):
        plt.plot(wavelengths, HDAC2_FS[i], label = HDAC2_label[i])
    
    
    plt.legend()
    plt.xlabel('Wavelength / nm')
    plt.ylabel('Lifetime / ns')
    plt.xlim([400, 700])
    plt.ylim([0.5, 2.0])
    
    plt.savefig(output_directory+'/Overlay plots', bbox_inches='tight')
    plt.show()
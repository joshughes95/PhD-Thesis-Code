# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:08:34 2022

@author: tscw79
"""

from __future__ import division
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askdirectory
import scipy
from scipy.optimize import minimize
import statsmodels.api as sm
plt.rcParams.update({'font.size': 12})

######################## Input Variables #####################################

Absorbance_wavelength = 390

Reference = ['LightOx17 in Toluene']

Reference_QY = [0.67]

Reference_RI = [1.497]

Acetone = 1.3588
CCl4 = 1.4607
CHCl3 = 1.4459
DCM = 1.4244
DMSO = 1.479
EtOAc = 1.3720
EtOH = 1.3611
IPA = 1.3776
MeCn = 1.344
MeCy = 1.441
Toluene = 1.497

Sample_RI = DMSO
######################### Absorbance ##########################################
Results_directory = askdirectory(title='Select Results Directory')
Absorbance_directory = askdirectory(title='Select Absorbance Directory')

#empty DF lists
df_list_abs = []
sheet_names = []

# Loop through folders
for folder in os.listdir(Absorbance_directory):
    
    sheet_names.append(folder)

    df_abs = []
    df_conc = []
    plt.figure()
    i = 0                 # use this to cycle through labels in the legend
    
    # loop through files within folders
    for filename in os.listdir(Absorbance_directory+'\\'+folder):
        df_output = pd.DataFrame(columns=('Concentration', 'Abs')) # df for each conc
        if os.path.splitext(filename)[1] == '.csv': # ignores aligent files
            if os.path.splitext(filename)[0] != 'Results':
                if os.path.splitext(filename)[0] != 'Blank':
                    df_conc.append(os.path.splitext(filename)[0])
                    
                    # Read input data (600 data points)
                    df = pd.read_csv(Absorbance_directory+'\\'+folder+'\\'+filename)        
                
                    wavelengths = df.iloc[1:601,0].astype(float)
                    absorbance = df.iloc[1:601,1].astype(float)
                    
                    # correct for any discrepancy in the baseline
                    if Reference[0] in filename:
                        correction_factor = np.average(absorbance[:200])
                    else:
                        correction_factor = np.average(absorbance[:200])
                    
                    # plot for each conc
                    plt.plot(wavelengths[100:500], absorbance[100:500] - correction_factor, label=df_conc[i])
                    i += 1
                    
                    # find absorbance at 390 nm
                    location_390 = np.where(np.round(wavelengths)==390)
                    absorbance_390 = absorbance[location_390[0]] - correction_factor
                    print(absorbance_390.iloc[0])
                    df_abs.append(absorbance_390.iloc[0])
     
    plt.ylabel('Abs')
    plt.xlabel('Wavelength/nm')
    #plt.xlim(300, 700)
    plt.legend()
    plt.ylim(-0.01, 0.12)
    plt.savefig(Absorbance_directory+'\\'+folder+'\\Absorption_Overlayed.png')
    plt.show()
    
    df_output['Concentration'] = df_conc
    df_output['Abs'] = df_abs
    df_list_abs.append(df_output)
        
Excelwriter = pd.ExcelWriter(Results_directory+'\\Absorbance.xlsx',engine="xlsxwriter")
for i, df in enumerate(df_list_abs):
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)

Excelwriter.save()


############################ Fluorometer ######################################

Fluorescence_directory = askdirectory(title='Select Fluorescence Directory')

def trapezium_rule(x, y):
    x = np.array(x)
    y = np.array(y)
    
    area = 0
    for j in range(len(x)-1):
        z = (x[j+1]-x[j])*(y[j+1]+y[j])/2
        if z > 0:
            area += z
        else:
            area = area
        #area += z
    return area

df_list_fluo = []
sheet_names = []

for folder in os.listdir(Fluorescence_directory):
    sheet_names.append(folder)
    for filename in os.listdir(Fluorescence_directory+'\\'+folder):
        if os.path.splitext(filename)[1] == '.csv':
            if '3D' in os.path.splitext(filename)[0]:
                df = pd.read_csv(Fluorescence_directory+'\\'+folder+'\\'+filename, skiprows=1)    
                df_conc = df.columns[2:]
                
                df_output = pd.DataFrame(columns=('Concentration', 'Fluorescence'))
                df_output['Concentration'] = df_conc
                df_fluo = []
            
                wavelengths = df[df.columns[0]][1:402].astype(float)
                Blank = df[df.columns[1]][1:402].astype(float)
                Blank_Fluorescence = trapezium_rule(wavelengths, Blank)
                plt.figure()
                
                for conc in df_conc:
                    emission = df[conc][1:402].astype(float)
                    fluorescence = trapezium_rule(wavelengths, emission)
                    fluorescence = fluorescence - Blank_Fluorescence
                    plt.plot(wavelengths, emission - Blank, label=conc)
                    print(fluorescence)
                    df_fluo.append(fluorescence)
                
                df_output['Fluorescence'] = df_fluo
                df_list_fluo.append(df_output)
                
                plt.xlabel('wavelength/nm')
                plt.ylabel('Fluorescence')
                plt.legend(loc='center right', bbox_to_anchor=(1.65, 0.5), 
                           ncol=1, fancybox=True, shadow=False)
                plt.savefig(Fluorescence_directory+'\\'+folder+'\\Fluorescence_Overlayed.png', bbox_inches='tight')
                plt.show()
        
Excelwriter = pd.ExcelWriter(Results_directory+'\\Fluorescence.xlsx',engine="xlsxwriter")
for i, df in enumerate(df_list_fluo):
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)

Excelwriter.save()

############################### Linear Regression #############################

cols = ['Compound', 'Gradient', 'Grad_err', 'Intercept', 'Int_err', 'R_sq']

Compounds = []
Gradients = []
Grad_errs = []
Intercepts = []
Int_errs = []
R_sqs = []

plt.figure()
colours = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'black', 'yellow', 'gray', 'pink']
for i, (absorbance, fluor) in enumerate(zip(df_list_abs,df_list_fluo)):  

    xdata = 1-10**(-absorbance['Abs'])
    ydata = fluor['Fluorescence']

    ols = sm.OLS(ydata, sm.add_constant(xdata))
    ols_result = ols.fit()
    coeff = ols_result.params
    std_err = ols_result.bse
    r_squared = ols_result.rsquared
    
    Compounds.append(sheet_names[i])
    Gradients.append(coeff[1])
    Grad_errs.append(std_err[1])
    Intercepts.append(coeff[0])
    Int_errs.append(std_err[0])
    R_sqs.append(r_squared)
    
    plt.plot(xdata, ydata, linestyle='', marker='x', label=sheet_names[i], color=colours[i])
    plt.plot(xdata, xdata*coeff[1]+coeff[0], color=colours[i])

plt.xlabel('Abs')
#plt.ylim(-5, 500)
plt.ylabel('Integrated Fluorescence')
plt.legend(loc='center right', bbox_to_anchor=(1.55, 0.5), 
           ncol=1, fancybox=True, shadow=False)
plt.savefig(Results_directory+'\\abs_v_Fluor.png', bbox_inches='tight')
plt.show()
    
df_output = pd.DataFrame(data={cols[0]: Compounds, cols[1]: Gradients, 
                               cols[2]: Grad_errs, cols[3]: Intercepts, cols[4]: Int_errs,
                               cols[5]: R_sqs})

Excelwriter = pd.ExcelWriter(Results_directory+'\\Gradients.xlsx',engine="xlsxwriter")
df_output.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
Excelwriter.save()

############################ Quantum Yields Calculations ######################
for j, reference in enumerate(Reference):
    def Quantum_yield(Ref, known_grad, unknown_grad, reference_RI, sample_RI):
        return Ref*(unknown_grad/known_grad)*((sample_RI/reference_RI)**2)
    
    def Quantum_yield_err(Ref, known_grad, unknown_grad, known_grad_err, unknown_grad_err, reference_RI, sample_RI):
        value = Quantum_yield(Ref, known_grad, unknown_grad, reference_RI, sample_RI)
        known_err = value - Quantum_yield(Ref, known_grad + known_grad_err, unknown_grad, reference_RI, sample_RI)
        unknown_err = value - Quantum_yield(Ref, known_grad, unknown_grad + unknown_grad_err, reference_RI, sample_RI)
        return np.sqrt(known_err**2 + unknown_err**2)
    
    Compounds = []
    QYs = []
    QYs_err = []
    QYs_perc_err = []
    
    ref_grad = df_output[df_output['Compound'] == reference]
    index = ref_grad.index.values[0]
    
    for i, comp in enumerate(df_output['Compound']):
        comp_grad = df_output[df_output['Compound'] == comp]
        if comp_grad.loc[i,'Compound'] != reference:
            print(ref_grad.loc[index,'Gradient'])
            comp_QY = Quantum_yield(Reference_QY[j], ref_grad.loc[index,'Gradient'], 
                                    comp_grad.loc[i,'Gradient'], Reference_RI[j], Sample_RI)
            comp_QY_err = Quantum_yield_err(Reference_QY[j], ref_grad.loc[index,'Gradient'], 
                                            comp_grad.loc[i,'Gradient'], ref_grad.loc[index,'Grad_err'], 
                                            comp_grad.loc[i,'Grad_err'], Reference_RI[j], Sample_RI)
            Compounds.append(comp)
            QYs.append(comp_QY)
            QYs_err.append(comp_QY_err)
            QYs_perc_err.append((comp_QY_err/comp_QY)*100)
        
        else:
            Compounds.append(comp)
            QYs.append(Reference_QY[j])
            QYs_err.append(0)
            QYs_perc_err.append(0)
    
    df_QY = pd.DataFrame(columns=('Compound', 'Quantum Yield', 'Quantum Yield Error', 'Quantum Yield Percentage Error'))
    df_QY['Compound'] = Compounds
    df_QY['Quantum Yield'] = QYs
    df_QY['Quantum Yield Error'] = QYs_err
    df_QY['Quantum Yield Percentage Error'] = QYs_perc_err
    
    Excelwriter = pd.ExcelWriter(Results_directory+'\\Quantum Yield with '+reference+' reference.xlsx',engine="xlsxwriter")
    df_QY.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
    Excelwriter.save()
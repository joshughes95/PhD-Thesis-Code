# -*- coding: utf-8 -*-
# Import relevant modules
#pylab inline
from __future__ import division
import numpy
import csv
import os
import matplotlib.pyplot as plt
from skimage.measure import block_reduce
from matplotlib.animation import FuncAnimation, PillowWriter
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm

################### VARIABLES #################################################

input_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2021-01-23 JC-1\BOTH\Processed_profiles_red'
output_directory = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2021-01-23 JC-1\BOTH\Figures_profiles_red'
output_bullseye = r'C:\Users\tscw79\OneDrive - LightOx\PhD\FRAP\Experiments\Data\2021-01-23 JC-1\BOTH\Bullseye_red_adjusted'

pixels = 1064**2
ROI_radius = 20
number_of_ROIs = 4
number_of_repeats = len(os.listdir(input_directory))/number_of_ROIs
pre_bleach_frames = 5
Laser_powers = numpy.array([0,1,10,50])
n = 6

ROI_radius_um = 2.5
conversion_factor = ROI_radius_um / ROI_radius

###################### MOVING AVERAGE FUNCTION ###############################

def moving_avg(x, n):
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0)) 
    return (cumsum[n:] - cumsum[:-n]) / float(n)

####################### EXTRACT DATA FROM EXCEL ###############################

#Loop over each repeat
processed_data = []
for filename in os.listdir(input_directory):
    data = csv.reader(open(input_directory+'\\'+filename, 'r'), delimiter=",", quotechar='|')

    #Read CSV, convert to an array and transpose
    data_extracter = []
    for row in data:
        data_extracter.append(row)
    data_array = numpy.zeros((len(data_extracter)-1,len(data_extracter[0])))
    for i in range(len(data_extracter)-1):
        data_array[i] = data_extracter[i+1]
    
    data_sorted = numpy.transpose(data_array)
    
    time = data_sorted[0]
    
    #To obtain Intenisty at each radii
    #Scale each ROI with area (r^2)
    #Subtract previous reading ie intenisty @ r=6: intensity reading of (I(r=6) *6^2) - (I(r=5)*5^2)
    data_scaled = numpy.zeros((len(data_sorted)-1,len(data_sorted[0])))
    for i in range(len(data_sorted)-1):
        data_scaled[i] = data_sorted[i+1]*((i+1)**2)
    
    data_at_radius = numpy.zeros((len(data_scaled), len(data_scaled[0])))
    for i in range(len(data_scaled)):
        if i == 0:
            data_at_radius[i] = data_scaled[i]
        else:
            data_at_radius[-i] = (data_scaled[-i]-data_scaled[-(i+1)])/(((len(data_scaled)-i+1)**2)-((len(data_scaled)-i)**2))
 
######################### NORMALISE DATA ######################################

#Try to normal for each radii, if the average over 5 point is zero, normalise using the total average over all radii
    #Normalise data over first 5 data points, across all ROIs by scale for number of pixels (area size)
    sum_norm_intensity = 0
    total_area = 0
    for i in range(len(data_at_radius)):
        sum_norm_intensity_2 = 0
        total_area += 5*(((i+1)**2)-(i**2))
        for j in range(int(pre_bleach_frames)):
            sum_norm_intensity_2 += data_at_radius[i][j]*(((i+1)**2)-(i**2))
        sum_norm_intensity += sum_norm_intensity_2
    
    #Normalise data over first 5 data points
    data_normalised = numpy.zeros((len(data_at_radius),len(data_at_radius[0])))
    for i in range(len(data_at_radius)):
        if numpy.mean(data_at_radius[i,:5]) == 0:
            data_normalised[i] = data_at_radius[i]*total_area / sum_norm_intensity
        else:
            data_normalised[i] = data_at_radius[i] / numpy.mean(data_at_radius[i,:5])

    processed_data.append(data_normalised)

#Convert data to an array
data_processed_ = numpy.array(processed_data)

#Find minimum array size
array_sizes = numpy.zeros((len(processed_data)))
for m in range(len(processed_data)):
    array_sizes[m] = len(data_processed_[m])
min_array_size = numpy.min(array_sizes)

#Reduce array size to minimum in data set
data_processed = numpy.zeros((len(data_processed_),int(min_array_size), len(data_processed_[0][0])))
for m in range(len(data_processed_)):
    for i in range(int(min_array_size)):
        data_processed[m][i] = data_processed_[m][i]
        
#Average data
data_processed_average = numpy.zeros((number_of_ROIs,len(data_processed[0]),len(data_processed[0][0])))
data_processed_error = numpy.zeros((number_of_ROIs,len(data_processed[0]),len(data_processed[0][0])))

for m in range(int(number_of_ROIs)):
    data_to_be_averaged = data_processed[m::number_of_ROIs]
    for i in range(len(data_processed[0])):
        for j in range(len(data_processed[0][0])):
            data_processed_average[m][i][j] = numpy.mean(data_to_be_averaged[:,i,j])
            data_processed_error[m][i][j] = numpy.std(data_to_be_averaged[:,i,j])
            
################## SAVE DATA IN PANDAS ########################################
df_list = []
df_list_std = []
sheet_names = ['0%', '1%', '10%', '50%']
columnA = ['Time']
columnB = (str(numpy.linspace(1,80,80))).split()[1:81]
columns = columnA + columnB

for i in range(len(Laser_powers)):
    
    df_average = pd.DataFrame(columns=(columns))
    df_average['Time'] = time
    df_std = pd.DataFrame(columns=(columns))
    df_std['Time'] = time
    
    for j, column in enumerate(columnB):
        df_average[column] = data_processed_average[i][j]
        df_std[column] = data_processed_error[i][j]
    
    df_list.append(df_average)
    df_list_std.append(df_std)
    
Excelwriter = pd.ExcelWriter(output_directory+'\\Bullseye.xlsx',engine="xlsxwriter")
Excelwriter_std = pd.ExcelWriter(output_directory+'\\Bullseye_std.xlsx',engine="xlsxwriter")
    
for i, (df, df_std) in enumerate(zip(df_list, df_list_std)):
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)
    df_std.to_excel(Excelwriter_std, sheet_name=sheet_names[i], index=False)

Excelwriter.save()
Excelwriter_std.save()

################### PLOT BULSEYE DATA ########################################

for m in range(int(number_of_ROIs)):
    radius = numpy.arange(1,len(data_processed_average[0])+1,1)
    x = moving_avg(radius[1:], n)
    data_to_plot = numpy.transpose(data_processed_average[m])
    ys = block_reduce(data_to_plot, block_size=(2,1), func=numpy.mean, cval=numpy.mean(data_to_plot))
    data_to_plot_error = numpy.transpose(data_processed_error[m])
    fig = plt.figure()
    ax = plt.axes(xlim=(0, len(data_to_plot[0])), ylim=(0.0, 2.0))
    line, = ax.plot([], [], lw=3)
    plt.vlines(ymin=0, ymax=3, x=ROI_radius, color='indianred')
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.title('Laser Power = %d' %Laser_powers[m])
    plt.xlabel('Radius/pixels')
    plt.ylabel('Normalised Intensity/au')
    
    #plt.legend('Time = {}s'.%{time[0]})
    #plt.legend.set_animated(True)
    
    def init():
        line.set_data([], [])
        return line,
    
    def animate(i):
        y = moving_avg(ys[i][1:], n)
        line.set_data(x, y)
        line.set_label('Time = %d s' %time[2*i])
        plt.legend()
        return line,
    
    anim = FuncAnimation(fig, animate, init_func=init,
    
                                    frames=50, interval=20, blit=True)
    
    anim.save(output_directory+'\\'+'Laser_power_'+str(Laser_powers[m])+'.gif', writer=PillowWriter(fps=2))
    plt.show()
    
    
    plt.plot(x, moving_avg(ys[0][1:], n), label='Pre-photoactivation', lw=2.5)
    plt.plot(x, moving_avg(ys[3][1:], n), label='Post-photoactivation', lw=2.5)
    #plt.plot(x, moving_avg(ys[45][1:], n), label='100 secs post-photoactivation', lw=2.5)
    plt.vlines(ymin=0, ymax=3, x=ROI_radius, color='indianred')
    plt.hlines(y=1, xmin=0, xmax=time[-1], color='indianred')
    plt.ylim(0.0,2.0)
    plt.xlim(0,80)
    plt.title('Laser Power = %d' %Laser_powers[m])
    plt.xlabel('Radius/pixels')
    plt.ylabel('Normalised Intensity/au')
    plt.legend(loc=4)
    #plt.savefig(output_directory+'\\'+'Overlay'+str(Laser_powers[m])+'.png', bbox_inches='tight')
    plt.show()
    
############################ DIFFUSION MODELLING ##############################

############################ 1. FIT TO GAUSSIAN ###############################

# Define the Gaussian function
def Gauss(x, A, B):
    y = A*numpy.exp(-(x**2)/(2*B))+1
    return y

# Pandas array to hold Gaussian data
df_list_gaussian = []
columns_gaussian = ['time', 'fit_A', 'fit_A_error', 'fit_B', 'fit_B_error', 'R^2', 'sigma', 'sigma_error']

xdata = moving_avg(radius[1:], n)*conversion_factor

for i, df in enumerate(df_list):
    fit_As = []
    fit_Bs = []
    fit_Cs = []
    r_sqs = []
    sigmas = []
    
    fit_A_error = []
    fit_B_error = []
    fit_C_error = []
    sigma_error = []
    
    ydata_t0 = 0
    # ydata_t0 = 1 - moving_avg(numpy.array(df_list[i].iloc[pre_bleach_frames, 1:])[1:], n) # adjust for bleaching
    
    for t in range(len(time) - pre_bleach_frames):
        # ydata = moving_avg(numpy.array(df_list[3].iloc[t + pre_bleach_frames, 1:])[1:], n) # Smooth data
        ydata = moving_avg(numpy.array(df_list[i].iloc[t + pre_bleach_frames, 1:])[1:], n) + ydata_t0 # Smooth data
        
        boundary = ([-numpy.inf, 0], [numpy.inf, numpy.inf])
        parameters, covariance = curve_fit(Gauss, xdata, ydata, p0=[-0.5,30], bounds=boundary)
        y_pred = Gauss(xdata, *parameters)
        r_sq = r2_score(ydata, y_pred)
        errors = numpy.sqrt(numpy.diag(covariance))
  
        fit_A = parameters[0]
        fit_B = parameters[1]
        #fit_C = parameters[2]
        
        sigma = numpy.sqrt(fit_B)
        
        fit_As.append(fit_A)
        fit_Bs.append(fit_B)
        r_sqs.append(r_sq)
        sigmas.append(sigma)
        #fit_Cs.append(fit_C)
        
        fit_A_error.append(errors[0])
        fit_B_error.append(errors[1])
        sigma_error.append(sigma - numpy.sqrt(fit_B - errors[1]))
        #fit_C_error.append(errors[2])
    
    df_gaussian = pd.DataFrame(data={columns_gaussian[0]: time[pre_bleach_frames:], columns_gaussian[1]: fit_As,
                                     columns_gaussian[2]: fit_A_error, columns_gaussian[3]: fit_Bs,
                                     columns_gaussian[4]: fit_B_error, columns_gaussian[5]: r_sqs, 
                                     columns_gaussian[6]: sigmas, columns_gaussian[7]: sigma_error})
    
    df_list_gaussian.append(df_gaussian)
    
Excelwriter = pd.ExcelWriter(output_bullseye+'\\Gaussian.xlsx',engine="xlsxwriter")
    
for i, df in enumerate(df_list_gaussian):
    df.to_excel(Excelwriter, sheet_name=sheet_names[i], index=False)

Excelwriter.save()

############################ 2. PLOT GAUSSIAN ################################
# Variables to change
time_points_to_plot = [0, 9, 17, 26, 34]
colours = ['blue', 'orange', 'green', 'purple', 'black']


for i, (df_G, df) in enumerate(zip(df_list_gaussian, df_list)):
    plt.figure()
    for j, t in enumerate(time_points_to_plot):
        
        fit_y = Gauss(xdata, df_G['fit_A'][t], df_G['fit_B'][t])
        plt.plot(xdata, moving_avg(numpy.array(df.iloc[t+pre_bleach_frames, 1:])[1:], n), color=colours[j], alpha=0.5)
        plt.plot(xdata, fit_y, '-', label='t={:.0f}s post-irradiation'.format(round(time[t],0)), color=colours[j])
    
    plt.vlines(ymin=0, ymax=3, x=ROI_radius*conversion_factor, color='indianred')
    plt.legend()
    plt.ylim(0.0,2.0)                           # match to plotted data above
    plt.xlabel('Radius / 'r'$\mu$m')
    plt.ylabel('Normalised Intensity / au')
    plt.savefig(output_bullseye+'\\Gaussian'+str(Laser_powers[i])+'.png', bbox_inches='tight')
    plt.show()

############################ 3. MODEL w^2 v t ################################

columns_D = ['Laser Power','Gradient', 'Gradient_error', 'Intercept', 'Intercept_error', 'Rsq', 'D', 'D_error']

tdata_1 = time[pre_bleach_frames+2:31].reshape((-1, 1))
tdata_2 = time[30:56].reshape((-1, 1))

Gradients_1 = []
Intercepts_1 = []
r_sqs_1 = []
Ds_1 = []

Gradients_error_1 = []
Intercepts_error_1 = []
Ds_error_1 = []

Gradients_2 = []
Intercepts_2 = []
r_sqs_2 = []
Ds_2 = []

Gradients_error_2 = []
Intercepts_error_2 = []
Ds_error_2 = []

for i, df in enumerate(df_list_gaussian):    

    ydata_1 = df['fit_B'][2:26]
    ydata_err_1 = df['fit_B_error'][2:26]

    # Run model
    model_1 = LinearRegression().fit(tdata_1, ydata_1)
    gradient_1 = model_1.coef_
    intercept_1 = model_1.intercept_
    r_sq_1 = model_1.score(tdata_1, ydata_1)
    D_1 = 0.5 * gradient_1
    
    ols_1 = sm.OLS(ydata_1, sm.add_constant(tdata_1))
    ols_result_1 = ols_1.fit()
    coeff_1 = ols_result_1.params
    std_err_1 = ols_result_1.bse
    r_squared_1 = ols_result_1.rsquared
    
    
    # DataFrame Entry
    Gradients_1.append(gradient_1)
    Intercepts_1.append(intercept_1)
    r_sqs_1.append(r_sq_1)
    Ds_1.append(D_1)
    
    Gradients_error_1.append(std_err_1['x1'])
    Intercepts_error_1.append(std_err_1['const'])
    Ds_error_1.append(0.5 * std_err_1['x1'])
    
    ydata_2 = df['fit_B'][25:51]
    ydata_err_2 = df['fit_B_error'][25:51]

    # Run model
    model_2 = LinearRegression().fit(tdata_2, ydata_2)
    gradient_2 = model_2.coef_
    intercept_2 = model_2.intercept_
    r_sq_2 = model_2.score(tdata_2, ydata_2)
    D_2 = 0.5 * gradient_2
    
    ols_2 = sm.OLS(ydata_2, sm.add_constant(tdata_2))
    ols_result_2 = ols_2.fit()
    coeff_2 = ols_result_2.params
    std_err_2 = ols_result_2.bse
    r_squared_2 = ols_result_2.rsquared
    
    
    # DataFrame Entry
    Gradients_2.append(gradient_2)
    Intercepts_2.append(intercept_2)
    r_sqs_2.append(r_sq_2)
    Ds_2.append(D_2)
    
    Gradients_error_2.append(std_err_2['x1'])
    Intercepts_error_2.append(std_err_2['const'])
    Ds_error_2.append(0.5 * std_err_2['x1'])
    
    # Plot 
    fit_y_1 = tdata_1*gradient_1 + intercept_1
    fit_y_2 = tdata_2*gradient_2 + intercept_2
    plt.errorbar(x=time[pre_bleach_frames+2:56]-time[pre_bleach_frames], y=df['fit_B'][2:51], yerr=df['fit_B_error'][2:51], marker='o', linestyle='', label='data', capsize=3)
    plt.plot(tdata_1-time[pre_bleach_frames], fit_y_1, '-', label='fit 1')
    plt.plot(tdata_2-time[pre_bleach_frames], fit_y_2, '-', label='fit 2')
    plt.xlabel('time / s')
    plt.ylabel(r'$\omega$$^2$ / 'r'$\mu$m$^2$')
    plt.annotate('Gradient = {} \u00B1 {}'.format(round(coeff_1['x1'],3), round(std_err_1['x1'],3)), xy=(0.1,0.9),xycoords='axes fraction')
    plt.annotate('R$^2$ = {:.3f}'.format(round(r_squared_1,3)), xy=(0.1,0.8),xycoords='axes fraction')
    plt.annotate('Gradient = {} \u00B1 {}'.format(round(coeff_2['x1'],3), round(std_err_2['x1'],3)), xy=(0.1,0.7),xycoords='axes fraction')
    plt.annotate('R$^2$ = {:.3f}'.format(round(r_squared_2,3)), xy=(0.1,0.6),xycoords='axes fraction')  
    plt.savefig(output_bullseye+'\\Diffusion'+str(Laser_powers[i])+'.png', bbox_inches='tight')
    plt.show()
    
df_diffusion = pd.DataFrame(data={columns_D[0]: sheet_names, columns_D[1]: Gradients,
                                 columns_D[2]: Gradients_error, columns_D[3]:Intercepts,
                                 columns_D[4]: Intercepts_error, columns_D[5]: r_sqs, 
                                 columns_D[6]: Ds, columns_D[7]: Ds_error})

Excelwriter = pd.ExcelWriter(output_bullseye+'\\Diffusion.xlsx',engine="xlsxwriter")
    
df_diffusion.to_excel(Excelwriter, index=False)

Excelwriter.save()

############################ 4. PLOT log(A) v log(t + t0) #################################


columns_D = ['Laser Power','Gradient', 'Gradient_error', 'Intercept', 'Intercept_error', 'Rsq', 'd', 'd_error']

Gradients = []
Intercepts = []
r_sqs = []
ds = []

Gradients_error = []
Intercepts_error = []
ds_error = []

for i, df in enumerate(df_list_gaussian):    

    ydata = df_list_gaussian[3]['fit_A'][2:50]
    for j, val in enumerate(ydata):
        if val < 0:
            ydata[j] = -ydata[j]
    ydata_log = numpy.log(ydata)
    ydata_err = df_list_gaussian[3]['fit_A_error'][2:50]
    ydata_err_log = abs(ydata_log - numpy.log(ydata + ydata_err))
    
    t0 = df_diffusion.iloc[3]['Intercept'] / (2 * df_diffusion.iloc[3]['D'])
    tdata = numpy.log((time[pre_bleach_frames+2:55] -time[pre_bleach_frames] + t0).reshape((-1, 1)))
    tdata = tdata[:]
    # Run model
    model = LinearRegression().fit(tdata, ydata_log)
    gradient = model.coef_
    intercept = model.intercept_
    r_sq = model.score(tdata, ydata_log)
    d = -2 * gradient
    
    ols = sm.OLS(ydata, sm.add_constant(tdata))
    ols_result = ols.fit()
    coeff = ols_result.params
    std_err = ols_result.bse
    r_squared = ols_result.rsquared
    
    
    # DataFrame Entry
    Gradients.append(gradient)
    Intercepts.append(intercept)
    r_sqs.append(r_sq)
    ds.append(d)
    
    Gradients_error.append(std_err['x1'])
    Intercepts_error.append(std_err['const'])
    ds_error.append(2 * std_err['x1'])
    
    # Plot 
    fit_y = tdata*gradient + intercept
    plt.errorbar(x=tdata-numpy.log(time[pre_bleach_frames]), y=ydata_log, yerr=ydata_err_log, marker='o', linestyle='', label='data', capsize=3)
    plt.plot(tdata-numpy.log(time[pre_bleach_frames]), fit_y, '-', label='fit')
    plt.xlabel('log(t+$t_{0}$)')
    plt.ylabel('log A')
    plt.annotate('Gradient = {} \u00B1 {}'.format(round(coeff['x1'],2), round(std_err['x1'],2)), xy=(0.1,0.9),xycoords='axes fraction')
    plt.annotate('R$^2$ = {:.2f}'.format(round(r_squared,2)), xy=(0.1,0.8),xycoords='axes fraction')    
    plt.savefig(output_bullseye+'\\Dimensionality'+str(Laser_powers[2])+'.png', bbox_inches='tight')
    plt.show()
    
df_dimensionality = pd.DataFrame(data={columns_D[0]: sheet_names, columns_D[1]: Gradients,
                                 columns_D[2]: Gradients_error, columns_D[3]:Intercepts,
                                 columns_D[4]: Intercepts_error, columns_D[5]: r_sqs, 
                                 columns_D[6]: ds, columns_D[7]: ds_error})

Excelwriter = pd.ExcelWriter(output_bullseye+'\\Dimensionality.xlsx',engine="xlsxwriter")
    
df_dimensionality.to_excel(Excelwriter, index=False)

Excelwriter.save()
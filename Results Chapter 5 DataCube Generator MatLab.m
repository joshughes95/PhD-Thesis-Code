%%test analysis script
%% SETUP PERAMITERS
% select the 1st folder to be analysed - note windows order does not
% always match matlab! - try hitting "run sectition" first to load the
% folder names and open the variable "all_files" to see the matlab order!
% the folder selected in the popup should be one level above this folder.
startingFileNumber = 21;

%setnumberofrows/colums - to match the data as recorded, if looking at
%multiple files that are not tile data just set numberOfRows to equal the
%number of data folders to analyse
numberofRows = 1;
numberofColums = 1;

%Set name for output toplevel folder where data is saved
analysisType = '\Analysis'; %use sample name if doing one measurement at a time
                                                            
AlphaMask = 1; % set to 0 to plot with no Aplha masking applied - TO BE IMPLEMENTED
AlphaScalefactor = 2.5; % Scaling factor for alpha mask
oneMinusAlpha = 0; %set to 1 for a 1-Alpha plot

%scaling for Intesnity normalisation, increase if needed
scalingFactorIntensity = 1.2;

%set the wavelength range to look over (need to convet from wavelength 1 =
%500 nm, 512 = 780 nm)
firstSensorPixel = 1;  % convert to wavelength by 0.6471 * lowest wavelength * pixelnumber 
lastSensorPixel = 512;
numberOfwavelengthstoPlot = 5; % number of wavelengths to fit beteen 500 and 720nm, evenly spaced

% set number of spectral pixels for moving sum, increase if data noisy
spectral_pixel_span = 32; 

% set threashold for when lifetime is set to NaN, based on peak bin fitted
count_threshold = 200;
mask_threshold = count_threshold;

% if you want to only want to create the lifetime data cubes:
% numberOfwavelengthstoPlot = 512,firstSensorPixel = 1 ,lastSensorPixel =
% 512 then set all the ploting / video options below to 0

plotImages = 1; % set to 1 to plot and save lifetime images set to 0 to simply save data
plotNormalisedImages = 0; % set to 1 to plot and save normalised lifetime images set to 0 to simply save data
createVideo = 1; % set to 1 to create a video of the computed lifetime images with histograms
videoQuality = 60; % set between 1 and 100 if output video too large
frameRate=2; % 45-60 works well for full spectral

%select 1st and last bin for fitting
binToFit1 = 10;
binToFit2 = 13;

lifetimeLow = 0.5; % for stained, 0.7,  1.5 for unstained, MHorrick 1
lifetimeHigh = 3.0; % for stained, 1.7,  2.8 for unstained, MHorrricks 2
      
%%

% Load file path and find number of folders - 1 level deap to workspaces!!
filePath = uigetdir;

newFolderIntesnity = [filePath, analysisType, '\Intensity'];
mkdir(newFolderIntesnity);
newFolderIntesnityNormalised = [filePath, analysisType, '\Intensity\Normalised'];
mkdir(newFolderIntesnityNormalised);
newFolderHistograms = [filePath, analysisType, '\Histograms'];
mkdir(newFolderHistograms);
newFolderMeanTau = [filePath, analysisType, '\Histograms\meanTau'];
mkdir(newFolderMeanTau);
newFolderLifetimeData = [filePath, analysisType, '\Lifetime_Data'];
mkdir(newFolderLifetimeData);

pixeldivider = round((lastSensorPixel-firstSensorPixel)/numberOfwavelengthstoPlot);
lastSensorPixel = pixeldivider*numberOfwavelengthstoPlot;

startingWavelength = 383.44;
Wavelength = 0;
for i = firstSensorPixel:pixeldivider:lastSensorPixel
    Wavelength = Wavelength +1;
    Wave  = round(i*0.6471 + startingWavelength);
    
    if plotImages == 1
    newFolderLifetimeImages = [filePath, analysisType,  '\Lifetime_', num2str(Wave), 'nm'];
    mkdir(newFolderLifetimeImages);
    
    newFolderHistogramsData = [filePath, analysisType, '\Histograms\', num2str(Wave),'nm'];
    mkdir(newFolderHistogramsData);
    end
    
    if plotNormalisedImages == 1
    newFolderLifetimeImages = [filePath, analysisType,  '\LifetimeNormalised_', num2str(Wave), 'nm'];
    mkdir(newFolderLifetimeImages);
    end 


end

pause(0.1)

all_files = dir([filePath, '\HistMode_no_pixel_binning\']);
all_files = struct2table(all_files);
all_files = all_files.name;
all_files = natsortfiles(all_files);
all_files = all_files(3:end);
numberofFiles = size(all_files,1);
load 'cmap.mat';

%%
%loop through images in date order, assumes the data was recorded in row by
%row with the same starting point, data saved to the corresponding folder,
%data workspaces must be only 1 level deap from the main folder
allIntensityImages={};
lifetimeImageData={};
lifetimeAlphaData={};
bin = 2;
currentFolder = pwd;
row = 0;

for r = 1:numberofRows
    row = row + 1;

    colum = 0;
    for k = 1:numberofColums
        colum = colum + 1;
        disp('row')
        disp(row)
        disp('column')
        disp(colum)
       
        fileNumber = row+colum-1 + ((row-1)*(numberofColums-1))+(startingFileNumber-1);
        imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
        %move to and load worspace from 1st subfolder
        currDir = [filePath,'\HistMode_no_pixel_binning\', cell2mat(all_files(fileNumber))];
        cd(currDir)
        disp('Loading workspace for folder:')
        disp(cell2mat(all_files(fileNumber)))
        load('workspace.frame_1.mat')
        sample = cell2mat(all_files(fileNumber));
        %return to matlab scripts directory
        cd(currentFolder)
        %set hitmode
        histMode = HIST_MODE;
        
        disp('Producing Intensity Image')

        % Produce and save intensity images
        [intensity_image] = Intensity_Image_Summation(bins_array_3, frame_size_x);
        climit = 'auto';
        plotter(intensity_image, newFolderIntesnity, sample, climit)
        
        allIntensityImages{row+colum-1 + ((row-1)*(numberofColums-1))} = intensity_image;
               
        %Calculate wavelength axis
        [wavelengths,wavenumbers] = Wavelength_Calculator(startingWavelength);
        
        disp('Performing Lifetime Calculations')
        % do lifetime fit calculations
        [parameters_cpu, selected_data_for_subtraction, bins_array_movsum_selected_reshaped] = test_LM_fitting_linear_gw(bins_array_3, histMode, spectral_pixel_span, binToFit1, binToFit2, frame_size_x);
        
        if plotImages == 1
            disp('Producing Lifetime Images and Plots for wavelength:')
        else
            disp('Producing Datacubes')
        end
        % Produce lifetime plots and histograms for various wavelengths
        
        numberofbins = size(selected_data_for_subtraction(:,1),1);
        selected_data_for_subtractionPeakbin = selected_data_for_subtraction(numberofbins,:);
        mask = selected_data_for_subtractionPeakbin;
        mask(mask<count_threshold)=0;
        mask(mask>count_threshold)=1;
        parameters_cpu(2,:) = parameters_cpu(2,:).*mask;
        tauLeastSquaresCPU = parameters_cpu(2,:); 
        tauLeastSquaresCPU(tauLeastSquaresCPU>4)=0;
        tauLeastSquaresReshaped = reshape(tauLeastSquaresCPU, [512 frame_size_x frame_size_x]);
        AlphaDataAll = reshape(selected_data_for_subtractionPeakbin, [512 frame_size_x frame_size_x]);
        
        
        wavelengthnumber = 0;
        for i = firstSensorPixel:pixeldivider:lastSensorPixel
            
           
            wavelengthnumber = wavelengthnumber+1;
            spectral_pixel = i;
            
            if plotImages == 1
            Currentwavelength = i*0.6471 + startingWavelength;
            disp(Currentwavelength)
            Lifetime_Image_Creation(spectral_pixel, bin, mask_threshold, sample, count_threshold, selected_data_for_subtraction, parameters_cpu, bins_array_movsum_selected_reshaped, wavelengths, lifetimeLow, lifetimeHigh, filePath, row, colum,analysisType, AlphaScalefactor, oneMinusAlpha, cmap, startingWavelength, frame_size_x);
            end
                
            if plotImages == 0
            tauLeastSquaresReshapedDisplayFrame = reshape(tauLeastSquaresReshaped(spectral_pixel,:,:),[frame_size_x  frame_size_x ]);
            array_movsum_selected = reshape(bins_array_movsum_selected_reshaped, size(bins_array_movsum_selected_reshaped, 2), size(bins_array_movsum_selected_reshaped, 1), size(bins_array_movsum_selected_reshaped, 3));
            bins_array_alpha = reshape(array_movsum_selected(spectral_pixel, bin, :),[frame_size_x frame_size_x]);
            end

        end
            lifetimeImageData{imageNumber} = tauLeastSquaresReshaped;
            lifetimeAlphaData{imageNumber} = AlphaDataAll;

    end
           


end

%% plot normalized intensity image - can be run independently if "allIntensityImages"
% is loaded and and the parameter initialisation section of the script is
% run
if plotNormalisedImages == 1
    AllnormalisationValue =[];
    disp('Creating Normalised Intesnity Images');
    row = 0; 
    for l = 1:numberofRows
        row = row + 1;
        colum = 0;
        for k = 1:numberofColums
            colum = colum + 1;
            imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
            ImagetoData =allIntensityImages{imageNumber};
            normalisationValue= max(max(ImagetoData));
            AllnormalisationValue(imageNumber) = normalisationValue;
        end
    end
    overallNormalisationValue = max(AllnormalisationValue);
    overallNormalisationValue = overallNormalisationValue/scalingFactorIntensity; 
    row = 0; 
    for l = 1:numberofRows
        row = row + 1;
        colum = 0;
        for k = 1:numberofColums
            colum = colum + 1;
            imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
            ImagetoPlot =allIntensityImages{imageNumber};
            IntensityImagesNormalised=ImagetoPlot/overallNormalisationValue;
            AllIntensityImagesNormalised{imageNumber}= IntensityImagesNormalised;
            climit = [0 1];
            plotter(IntensityImagesNormalised, newFolderIntesnityNormalised, row, colum, climit);
        end
    end

    save([newFolderIntesnityNormalised,'\AllIntensityData.mat'],'allIntensityImages')
    save([newFolderIntesnityNormalised,'\AllIntensityImagesNormalised.mat'],'AllIntensityImagesNormalised')
    else
end
%% Plot global normalized Alpha Lifetime plots, - can be run independently 
% if "LifetimeAlphaData"  and "LifetimeImageData" are loaded  along 
% with runing the  initialiation section of the script

if plotNormalisedImages == 1
    disp('Creating Normalised Lifetime Images');
    NormalisedAlphaData=[];
    for z = firstSensorPixel:pixeldivider:lastSensorPixel
        row = 0;
        for l = 1:numberofRows
            row = row + 1;
            colum = 0;

            for k = 1:numberofColums
                colum = colum + 1;
                imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
                imageAlphaData = lifetimeAlphaData{imageNumber};
                AlphaDataWavelength =imageAlphaData(z,:,:);
                normalisationValue= max(max(AlphaDataWavelength));
                NormalisedAlphaData(z,imageNumber) = normalisationValue;
            end
        end
    end
    

        a=0;
    for z = firstSensorPixel:pixeldivider:lastSensorPixel
        row = 0;
        a=a+1;
        CurrentWavelength = round(z*0.6471 + startingWavelength);
        OverallNormalisedAlphaFactor = max(NormalisedAlphaData(z,:));
        for l = 1:numberofRows
            row = row + 1;
            colum = 0;

            for k = 1:numberofColums
                colum = colum + 1;
                imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
                ImageData =lifetimeImageData{imageNumber};
                AlphaData = lifetimeAlphaData{imageNumber};
                ImagetoPlot = squeeze(ImageData(z, : , :));

                if oneMinusAlpha == 1
                AlphaData = squeeze(AlphaData(z,:,:));
                AlphatoPlot = AlphaScalefactor*(1-(AlphaData/OverallNormalisedAlphaFactor));
                folder = [filePath,analysisType, '\LifetimeNormalised_', num2str(CurrentWavelength),'nm','\r', num2str(row),'_c', num2str(colum), '_Lifetime ', num2str(CurrentWavelength),'nm.tif'];
                lifetimeplotter(ImagetoPlot, folder, AlphatoPlot, lifetimeLow, lifetimeHigh, cmap)
                else
                AlphaData = squeeze(AlphaData(z,:,:));
                AlphatoPlot = AlphaScalefactor*squeeze(AlphaData)/OverallNormalisedAlphaFactor;
                folder = [filePath,analysisType, '\LifetimeNormalised_', num2str(CurrentWavelength),'nm','\r', num2str(row),'_c', num2str(colum), '_Lifetime ', num2str(CurrentWavelength),'nm.tif'];
                lifetimeplotter(ImagetoPlot, folder, AlphatoPlot, lifetimeLow, lifetimeHigh, cmap)
                end
            end
        end
    end
    else
end

%% Create Video 
% can run if "LifetimeAlphaData"  and "LifetimeImageData" are loaded  along 
% with runing the  initialiation section of the script

    
if createVideo ==1
    disp('Creating Video');
    load cmap
    % setup a video file to populate
    filename = fullfile(newFolderLifetimeData, '\processed_video.avi');
    compressedVideo = VideoWriter(filename); 
    %set the final video frame rate
    compressedVideo.FrameRate = frameRate;
    compressedVideo.Quality = videoQuality;
    % open the video file ready to start adding frames
    open(compressedVideo);

    NormalisedAlphaData=[];
    for z = firstSensorPixel:pixeldivider:lastSensorPixel
        row = 0;
        for l = 1:numberofRows
            row = row + 1;
            colum = 0;

            for k = 1:numberofColums
                colum = colum + 1;
                imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
                imageAlphaData = lifetimeAlphaData{imageNumber};
                AlphaDataWavelength =imageAlphaData(z,:,:);
                normalisationValue= max(max(AlphaDataWavelength));
                NormalisedAlphaData(z,imageNumber) = normalisationValue;
            end
        end
    end

        a=0;
    for z = firstSensorPixel:pixeldivider:lastSensorPixel
        row = 0;
        a=a+1;
        
        CurrentWavelength = round(z*0.6471 + startingWavelength);
        OverallNormalisedAlphaFactor = max(NormalisedAlphaData(z,:));
        for l = 1:numberofRows
            row = row + 1;
            colum = 0;

            for k = 1:numberofColums
                colum = colum + 1;
                imageNumber = row+colum-1 + ((row-1)*(numberofColums-1));
                ImageData =lifetimeImageData{imageNumber};
                ImagetoPlot =squeeze(ImageData(z,:,:));
                AlphaData = lifetimeAlphaData{imageNumber};

                if oneMinusAlpha == 1
                AlphaData = squeeze(AlphaData(z,:,:));
                AlphatoPlot = AlphaScalefactor*(1-(AlphaData/OverallNormalisedAlphaFactor));
                [compressedVideo] =  lifetimeplotter_video(ImagetoPlot,  AlphatoPlot, lifetimeLow, lifetimeHigh, compressedVideo,CurrentWavelength, bins_array_movsum_selected_reshaped, bin, mask_threshold, z, cmap);
                else
                AlphaData = squeeze(AlphaData(z,:,:));
                AlphatoPlot = AlphaScalefactor*AlphaData/OverallNormalisedAlphaFactor;
                [compressedVideo] =  lifetimeplotter_video(ImagetoPlot,  AlphatoPlot, lifetimeLow, lifetimeHigh, compressedVideo,CurrentWavelength, bins_array_movsum_selected_reshaped, bin, mask_threshold, z, cmap);
                end
            end
        end
    end
    %close the video file
    close(compressedVideo);
    else
end

%%
 % Impllement saving in non cell format?
% lifetimeImageDatatoSave = cell2mat(lifetimeImageData);
% lifetimeImageDatatoSave = permute(reshape(lifetimeImageData,[frame_size_x, numberOfwavelengthstoPlot, frame_size_x]),[2 1 3]);
% lifetimeAlphaDatatoSave = cell2mat(lifetimeAlphaData);
% lifetimeAlphaDatatoSave = permute(reshape(lifetimeAlphaData,[frame_size_x, numberOfwavelengthstoPlot, frame_size_x]),[2 1 3]);
disp('Saving Datacubes');
save([newFolderLifetimeData,'\LifetimeImageData.mat'],'lifetimeImageData')
save([newFolderLifetimeData,'\LifetimeAlphaData.mat'],'lifetimeAlphaData')

%create Meta Data File
xgalvo_step_size = 35;
firstWavelength = 0.6471 *firstSensorPixel + startingWavelength;
lastWavelength = 0.6471 *lastSensorPixel + startingWavelength;
AnalysedfolderName = split(filePath,"\");
metaDataFolderName = strcat('\metaData_', string(AnalysedfolderName(2)), '.csv');
metaData = {};
metaData{1}  = ["Folder Analysed " ,  AnalysedfolderName(2)];
metaData{2}  = ["Hist mode " ,  num2str(HIST_MODE)];
metaData{3}  = ["Pstop" , num2str(CODE_PSTOP)];
metaData{4}  = ["Frame Size " , num2str(frame_size_x)];
metaData{5}  = ["Step Size " , num2str(xgalvo_step_size)];
metaData{6}  = ["Exposure Time " , num2str(exposure_time_us)];
metaData{7}  = ["Bin1 for Fitting" , num2str(binToFit1)];
metaData{8}  = ["Bin2 for Fitting " , num2str(binToFit2)];
metaData{9}  = ["Alpha Mask Enabled " , num2str(AlphaMask)];
metaData{10} = ["Alpha Scale Factor " , num2str(AlphaScalefactor)];
metaData{11} = ["1-Alpha Enabled " , num2str(oneMinusAlpha)];
metaData{12} = ["Count Threashold " , num2str(count_threshold)];
metaData{13} = ["Number of Wavelengths Analysed " , num2str(numberOfwavelengthstoPlot)];
metaData{14} = ["Starting Wavelength " , num2str(firstWavelength)];
metaData{15} = ["Last Wavelength " , num2str(lastWavelength)];
metaData{16} = ["Short Lifetime for Plots " , num2str(lifetimeLow)];
metaData{17} = ["Long Lifetime for Plots " , num2str(lifetimeHigh)];
metaData{18} = ["Vidio Frame Rate " , num2str(frameRate)];
metaData{19} = ["Vidio Compression (%) " , num2str(videoQuality)];
metaData{19} = ["Moving Average Size " , num2str(spectral_pixel_span)];
metaData = splitvars(cell2table(metaData'));
metaData.Properties.VariableNames = {'Variable' 'Value'};
writetable(metaData,strcat(filePath, analysisType , metaDataFolderName));

pause (1);
clear all

%process_MMI_sample
%
% This script uses several custom functions to grab raw MMI files and plots

%%% THIS IS A SAMPLE SCRIPT!


%All Settings in this file are pre-set for the "111 - LearnToProcess" data
%% File Inputs
%%% sampleFile assumed to be Name_stuff..._FSn.tif (FSn for frequency ranges)
%%% ---'Name' will be used for tagging all files
sampleFile = 'ILAS0d1mm_120_All.tif';
FStag      = '_All';

%File variables must be full path
CalibFile  = 'Calibration_P1Under_P45W60.mat';
darkMFile  = 'Dark_120.tif';  %name of "Dark" image file. Usually only number is updated

samplename = strtok(sampleFile,'_');   %for plot and Save

%%% NOTE: Here, "samplename" simply grabs all the characters up to the
%%% first underscore. Here, it grabs "ILAS0d1mm" and names all output files with that name.


%%% These save and plot folders are never changed
saveFolder = '.\Processed\';
plotFolder = [saveFolder 'SamplePlots\' samplename '\'];


%% Processing Parameters

% Standard processing for a SET of data is to SAVE the first ROI info
% using "saveROI = 1" and don't load any pre-existing ROIs (loadROI = 0).
% Usually, the "saveROItag" will simply be the samplename.
%
% Then, after you process the first sample, set
% saveROI to 0 and loadROI to 1. The loadROItag should be the first
% sample's name (here, already done to 'ILAS0d1mm' for this example).
%
% That means... all your following samples will be analyzed with the SAME
% ROI, which is great for comparing data!

saveROItag     = samplename;
loadROItag     = 'ILAS0d1mm';
saveROI        = 1;          %BE MINDFUL OF THESE EACH TIME
loadROI        = 0;          %usually save OR load
%both save and load to ...\Processed\ROI_Data_roiTag.mat


%%% NOTE: The only parameter you're likely to use below is the "TroiS"
%%% parameter that sets the small ROI for 'trace' analysis. It's in Plot
%%% Parameters below.

% For now, you can ignore all the options below if you like!

%% Everything below this changes infrequently
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmlfile     = ['Samples' FStag '.xml'];  %parameter file name. Usually 'Samples_FSn.xml'

%% Make Folders
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end
if ~exist(plotFolder,'dir')
    mkdir(plotFolder);
end
%% Load Parameters
MeanRawReps = 1;  %set to 1 to average reps. not well coded for not-averaging reps
roiCrop     = 1;  
%1 to choose roi location and size,
%[w h] to set width and height but do location by hand,
%[left top w h] to set size and location programmatically
%% Plot Parameters
TroiS        = [25 25]'; 

%TroiS sets the ROI information when doing the trace analysis. "Trace"
%means an ROI is chosen for an image and then analyzed over all spatial
%frequencies.

% VALUES
%       1: User chooses size and location, but is always a square
%  [w, h]: User chooses location. Width 'w' and height 'h' are set. (in pixels)

% [Vstart,Hstart,Vlength,Hlength]:
% This format defines the entire ROI size and position.
%  Vstart:  defines the vertical start of the ROI
%  Hstart:  defines the horizontal start of the ROI
%  Vlength: defines the vertical length of the ROI
%  Hlength: defines the horizontal length of the ROI

%% RUN
%organize DEMOD info
demodGroup.demod = 1;
demodGroup.ACmean= 1;
demodGroup.acAdd3= 0;
demodGroup.ACROI = 0;
demodGroup.getOnePhase = 0;
demodGroup.AcqLoop = 0;
if strcmpi(darkMFile(end),'f')%.tif file, get data
    darkM = get_darkM(darkMFile,1);
else
    load(darkMFile,'darkM');  %.mat file, load data
end

%load parameter info
X = xml2struct(xmlfile);
numFreq = str2double(X.XML_Cluster.Projection.NumFreq.Text);
Freqs   = linspace(str2double(X.XML_Cluster.Projection.LowFreq.Text),...
          str2double(X.XML_Cluster.Projection.TopFreq.Text),...
          numFreq);
numPhi = str2double(X.XML_Cluster.Projection.NumPhi.Text);
numReps = str2double(X.XML_Cluster.Acquisition.RepInfo.NumReps.Text);

%load calibration and ROI data
load(CalibFile,'A','W');
if loadROI
    if strcmpi(loadROItag,''), roiStartTag = 'ROI_Data';
    else, roiStartTag = 'ROI_Data_'; end
    load(['.\Processed\' roiStartTag loadROItag '.mat']); roiCrop = roi; TroiS = Troi;
end
%Grab files... single for normal sample, three to demod for SI sample.

clearvars 'acmeanROI';
[Btemp,roi,badreps] = tif2Bmat(sampleFile,numReps,MeanRawReps,darkM,roiCrop,numFreq,numPhi,0);
B = zeros(size(Btemp,1),size(Btemp,2),4,4,numFreq);

ACmeanf   = 0;
count     = 0;
for f = 1:numFreq
    if Freqs(f)~=0  %if processing anything besides f == 0
        ACmeanf = 1; 
        if count == 0,  count = count +1;end
    end
    for Mrow = 1:4
        for Mcol = 1:4
            %ASSUMES 3 PHASES
            %                         if count == 1, acmeanROI = ACROI; count = 2; end
            [B(:,:,Mrow,Mcol,f),~] = miDemod(squeeze(Btemp(:,:,Mrow,Mcol,f,1)),...
                squeeze(Btemp(:,:,Mrow,Mcol,f,2)),...
                squeeze(Btemp(:,:,Mrow,Mcol,f,3)),ACmeanf,0);
        end
    end
end

clearvars Btemp;

%get M
M = b2m(B,A,W);

%% plot data
for F = 1:numFreq %Plotting Full Mueller Matrices of Images
    freqTag = [': Freq ' num2str(Freqs(F)) ' mm^-^1'];
    plotMimages(B(:,:,:,:,F));
    caxis auto;
    title(['B Matrix of ' samplename freqTag]);
    
    %Plot Options:  'crange',crange,'cmap',cmap,'fhandle',fhandle
    Fmmi = plotMimages(squeeze(M(:,:,:,:,F)),1);
    title(['Mueller Matrix of ' samplename freqTag ': normed']);
    saveas(Fmmi,[plotFolder samplename FStag '_' num2str(Freqs(F)) '_MMIs.fig']);
    saveas(Fmmi,[plotFolder samplename FStag '_' num2str(Freqs(F)) '_MMIs.png']);
end

Troi = zeros(4,size(TroiS,2));
Mmeans = zeros(4,4,length(Freqs),size(TroiS,2));
Mstds = zeros(4,4,length(Freqs),size(TroiS,2));
for R = 1:size(TroiS,2) %how many ROIs?
    [Troi(:,R), Mmeans(:,:,:,R), Mstds(:,:,:,R), Fsi] = plotMMIvsSI(M,Freqs,TroiS(:,R),1,0,1);
    
    if size(TroiS,2)>1, TitleTag = ...
            [' Polarization over spatial frequency - ' num2str(R)];
    else
        TitleTag = ' Polarization over spatial frequency';
    end
    title([samplename TitleTag]);
    saveas(Fsi,[plotFolder samplename FStag '_MMI_over_SI_' num2str(R) '.fig']);
    saveas(Fsi,[plotFolder samplename FStag '_MMI_over_SI_' num2str(R) '.png']);
    if R == size(TroiS,2)
        note = 'Mmeans is normalized by m11';
        save([saveFolder samplename FStag '_siTraceData.mat'],'Troi','TroiS','Mmeans','Mstds','note');
    end
end

%save
save([saveFolder samplename FStag],'B','M','A','W','roi','Freqs','Troi','MeanRawReps','darkM','badreps','demodGroup');
if saveROI, save([saveFolder 'ROI_Data_' saveROItag '.mat'],'Troi','roi'); end

% eval(['B' out_tag ' = B;']);
eval(['M_' samplename ' = M;']);
clearvars B M P out_tag phases tag ans Fmmi F MAll;


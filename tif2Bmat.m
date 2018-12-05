% [Bout,roi,badreps] = tif2Bmat(filename,reps(,meanFlag,darkM,ROICrop,numFreq,numPhi,AcqLoopOption,AvoidReps)
%
% This function will read a tif stack and return the measurement matrix B
% for every pixel. This is done in the context of Mueller Matrix Imaging.
%
% INPUT
% filename: filepath string to .tif file
% numReps:  number of repititions at each state
%
% OPTIONS -- If skipping, enter []
% meanFlag: (optional) set to 1 to average all reps. Assumes reps are taken
%           at each state (N state 1, then N of state 2, ...). Set to 0 to
%           output all reps (y,x,4,4,N). Default is 1.
%
% darkM:    (optional) give darkcount to subtract from every image (this is
%           assumed to be a scalar value). Default is no subtraction.
%
% ROICrop:  (optional) various ROI actions can be chosen ->
%           1x4 position vector: [left start, top start, width, height]
%           1x2 roi size vector: [width height] of UI placed ROI
%           1: flag for complete user freedom with ROI rect
%           Default is no cropping.
%
% numFreq:  (optional) For structured illumination, this is the number of
%           spatial frequencies to load in. It assumes the save order is
%           [rep, phase, freq, M state]. Default is 1.
%
% numPhi:   (optional) For structured illumination, this is the number of
%           spatial frequency phases to load in. It assumes the tiff order
%           is [rep, phase, freq, M state]. Default is 1. (usual is 3).
%
% AcqLoopOption: (optional) For structured illumination, are the Mueller
%           matrix states looped through while spatial frequency was held
%           constant, or vice versa?  0 (Default) means MM loops first. 1
%           is for spatial frequency (old data) loops first.
%
% AvoidReps: (optional) vector list of reps to avoid loading
%
% OUTPUT
% Bout:     Measurement matrix of images. This function assumes, like
%           get_Calib_data.m, a state sequence for acquisition. Each LC has
%           two 'positions' (0 and 1). Represent each state with the vector
%           L1L2L3L4, where L1 is the first LC the source goes through, and
%           L4 is the last. The sequence starts with 0000 and ends with
%           1111. The sequence mimics binary counting (0000,0001,0010...).
%           Matrix is organized as (row,col,Mrow,Mcol (,N))
% roi:      (optional) the roi that was provided or chosen to crop B0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Bout, varargout] = tif2Bmat(filename,numReps,varargin)
screenSize = get(0,'ScreenSize'); screenSize(1:2) = []; %width,height

%grab x,y dimensions of images. Assume the same for sample images.
iminfo = imfinfo(filename);
height = iminfo(1).Height;
width  = iminfo(1).Width;

numFreq = 1;
numPhi  = 1;
AcqLoopOption = 0;
AvoidReps = [];
if nargin > 5, numFreq = varargin{4}; end
if nargin > 6, numPhi = varargin{5}; end
if nargin > 7, AcqLoopOption = varargin{6}; end
if nargin > 8, AvoidReps = varargin{7}; end
%set defaults for optional inputs
meanFlag = 1;

%read optional arguments
if nargin > 2, meanFlag = varargin{1}; end;

%Read in data
if meanFlag, Bout = zeros(height,width,4,4,numFreq,numPhi); Itemp = zeros(height,width,numReps);
else         Bout = zeros(height,width,4,4,numFreq,numPhi,numReps); %#ok<SEPEX>
end
%check for darkM
if nargin > 3
    darkM = varargin{2};
else
    darkM = 0;
end
count = 0;
badreps = [];
sig_thresh = darkM+sqrt(darkM);
if AcqLoopOption  %Old Looping. Acq loop through SI first, then MM
    for col = 1:4
        for row = 1:4
            for f = 1:numFreq
                for p = 1:numPhi
                    for r = 1:numReps
                        %assumes save order changes rep -> phi -> freq -> m
                        count = count+1;
                        if meanFlag
                            Itemp(:,:,r) = imread(filename,count);
                            if Itemp(:,:,r) < sig_thresh, badreps = [badreps count]; end
                            if r == numReps
                                %Signal < Threshhold will not be considered
                                goodreps = mean(mean(Itemp(:,:,:),1),2)>sig_thresh;
                                
                                Bout(:,:,row,col,f,p) = mean(Itemp(:,:,goodreps),3);
                            end
                        else
                            Bout(:,:,row,col,f,p,r) = imread(filename,count);
                        end
                    end
                end
            end
        end
    end
else   %New looping. Acq loops through Mueller states first, then SI.
    for f = 1:numFreq
        for p = 1:numPhi
            for col = 1:4
                for row = 1:4                   
                    for r = 1:numReps
                        %assumes save order changes rep -> phi -> freq -> m
                        count = count+1;
                        if meanFlag
                            Itemp(:,:,r) = imread(filename,count);
                            if Itemp(:,:,r) < sig_thresh, badreps = [badreps count]; end
                            if r == numReps
                                %Signal < Threshhold will not be considered
                                goodreps = mean(mean(Itemp(:,:,:),1),2)>sig_thresh;
                                Bout(:,:,row,col,f,p) = mean(Itemp(:,:,goodreps),3);
                            end
                        else
                            Bout(:,:,row,col,f,p,r) = imread(filename,count);
                        end
                    end
                end
            end
        end
    end
end
%Subtract Dark value.  assigned above
if nargin > 3
    Bout = Bout - darkM;
end
if nargout == 3, varargout{2} = badreps; end
%ROI management
roi = [];
if nargin > 4
    if length(varargin{3}) == 4     %ROI provided
        roi = varargin{3};
    elseif length(varargin{3}) == 2 %ROI dimensions provided
        roiW = varargin{3}(1);
        roiH = varargin{3}(2);
        F = figure; imagesc(Bout(:,:,1,1,1,1,1));
        hr = imrect(gca,[width/2 height/2  roiW roiH]);
        setResizable(hr,false);
        hm = msgbox({['Drag the rectangle to choose the center of a '...
            num2str(roiW) 'x' num2str(roiH) ' pixel ROI, then'];...
            ['                                    ' ...
            'DOUBLE CLICK it to finish.']},'Help Me Help You');
        hm.Position = [screenSize(1)/3.2 screenSize(2)/2.8 290 60];
        uiwait(hm);
        roi = wait(hr); close(F);
    elseif varargin{3} == 1         %ROI must be square
        F = figure; imagesc(Bout(:,:,1,1,1,1,1));
        hr = imrect(gca,[width/3 height/3 width/4 height/4]);
        setFixedAspectRatioMode(hr,true);
        msgbox({'Move and resize the rectangle to choose your ROI, then';...
            ['                                    ' ...
            'DOUBLE CLICK it to finish.']},'Help Me Help You');
        roi = wait(hr); close(F);
    end
    %Crop it
    BoutC = imcrop(Bout(:,:,1,1,1,1,1),roi);
    BoutC = zeros(size(BoutC,1), size(BoutC,2),4,4,size(Bout,5),size(Bout,6),size(Bout,7));
    for f = 1:numFreq
        for p = 1:numPhi
            for r = size(Bout,7)
                for row = 1:4
                    for col = 1:4
                        BoutC(:,:,row,col,f,p,r) = imcrop(Bout(:,:,row,col,f,p,r),roi);
                    end 
                end
            end
        end
    end
    
    Bout = BoutC;
    if nargout >1, varargout{1} = roi; end
end

%% miDemod(I1,I2,I3)

%% This basic function does the phase stepping for standard 3-phase SFDI
%This function assumes each matrix has the first row and column as the
%image height and width. The function will iterate over each other
%dimension and output a matrix the same size as the input.

%Handles up to  4D  at the moment

%OPTION
% ACmean
%        this processing is for AC frequencies only. It take the mean of
% each image and sets them to be the same. This helps clean up the AC
% demodulation and does not bias the result in any way.

function [demod, varargout] = miDemod(I1, I2, I3, varargin)

if all(size(I1)==size(I2)) && all(size(I1)==size(I3))
    
    demod = zeros(size(I1));
    acmeanROI = [1 1 size(I1,2) size(I1,1)];

    if nargin > 3
        if varargin{1} == 1
            %handle ROI info
            if nargin > 4
                if length(varargin{2})>1
                    acmeanROI = varargin{2};
                elseif varargin{2} == 1
                    F = figure; imagesc(I1(:,:,1,1));
                    hr = imrect(gca,[20 20 20 20]);
                    msgbox({'Move and resize the rectangle to choose your ROI and';...
                        ['                                    ' ...
                        'AVOID glare. DOUBLE CLICK it to finish.']},'Help Me Help You');
                    acmeanROI = wait(hr); close(F);
                end
            end
            %do "ac mean" correction
            for i = 1:size(I1,3)
                for j = 1:size(I1,4)
                    mDCs = [mean(mean(imcrop(I1(:,:,i,j),acmeanROI),1),2) ...
                        mean(mean(imcrop(I2(:,:,i,j),acmeanROI),1),2) ...
                        mean(mean(imcrop(I3(:,:,i,j),acmeanROI),1),2)];
                    %take away current mean and add overall mean...
                    %this makes a constant DC offset between images
                    I1(:,:,i,j) = I1(:,:,i,j)+mean(mDCs)-mDCs(1);
                    I2(:,:,i,j) = I2(:,:,i,j)+mean(mDCs)-mDCs(2);
                    I3(:,:,i,j) = I3(:,:,i,j)+mean(mDCs)-mDCs(3);
                end
            end
        end
    end
    %demodulate
    for i = 1:size(I1,3)
        for j = 1:size(I1,4)
            demod(:,:,i,j) = sqrt(2)/3 * sqrt(...
                (I2(:,:,i,j)-I1(:,:,i,j)).^2 ...
                + (I3(:,:,i,j)-I2(:,:,i,j)).^2 ...
                + (I1(:,:,i,j)-I3(:,:,i,j)).^2);
        end
    end
    
else
    error('All three images must be the same size!');
end
if nargout>1, varargout{1} = acmeanROI; end
end
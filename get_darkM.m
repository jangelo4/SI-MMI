%% [darkM (darkMs)] = get_darkM(file(,plot_checks))
%
% This function is for MMI processing. Used in get_Calib_data.m and
% process_MMI_sample.m. It expects a 'dark' image stack (camera shutter
% does not open, but this is not necessary).  It averages over the entire
% image to a point. If a frame has an average value that is more than a
% standard deviation away from the other points, it is disregarded.
%
% INPUT: DarkFile - the dark tiff file
%        plot_checks - (optional) 1: plot the quality assurance checks to
%        see if they are done correctly. 0 turns this off. 0 is default.
%
% OUTPUT: darkM - the averaged dark value over all pixels and frames
%         darkMstds - standard deviation of dark images
%         darkMs - (optional) this is the value from each frame.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [darkM, varargout] = get_darkM(DarkFile, varargin)


iminfo = imfinfo(DarkFile);

reps_d = length(iminfo);
height = iminfo(1).Height;
width  = iminfo(1).Width;
%% Dark Prep
%read images
dark = zeros(height,width,reps_d);
for i = 1:reps_d
    dark(:,:,i) = imread(DarkFile,i);
end
darkMs = squeeze(mean(mean(dark,1),2));  %Find Average for each rep
darkM  = mean(darkMs);                   %Find Average for everything
goodDs = ones(reps_d,1);                 %assume everything is good!
%measure average of the standard deviation within each rep
darkstdM = mean(std(reshape(dark,height*width,reps_d)));
for i = 1:reps_d
    if darkMs(i)> darkstdM+darkM || darkMs(i)< darkM-darkstdM
        goodDs(i) = 0;
    end
end
%report list
if isempty(find(~goodDs,1))
    dBadList = 'none';
else
    dBadList = sprintf('%.0f, ' , find(~goodDs)); %find bad reps
    dBadList = dBadList(1:end-2);% strip final comma
end
hr = figure; plot(darkMs,'o-'); hold on; 
plot(ones(5,1)*(darkM-darkstdM),'k--'); %%show standard deviation lines
plot(ones(5,1)*(darkM+darkstdM),'k--');
ylim([darkM-1.25*darkstdM darkM+1.25*darkstdM]); %see within 1 standard deviation

title(['Averaged Dark Reps. Bad reps: ' dBadList]);
hm = msgbox({'Check the plot for dark values that are abnormal,';...
    '(greater than one standard deviation... ~sqrt(mean))';...
    'make sure they''re listed in then title, then click OK.'},'Find Acquisition Errors');
% hm.Position = [screenSize(1)/3.2 screenSize(2)/2.8 270 60];
uiwait(hm);
close(hr);
darkM = mean(mean(mean(dark(:,:,logical(goodDs)))));

if nargout == 2
    varargout{1} = darkstdM;
end
if nargout == 3
    varargout{2} = darkMs;
end


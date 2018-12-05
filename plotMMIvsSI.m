%[roi,MNmeans,MNstds f] = plotMMIvsSI(M,freqs (,roi, elements, errbars, rawm11, f) )
%
% This function plots elements of M vs spatial frequency
%
% INPUT: M - Mueller Matrix. Built (Irow,Icol,Mrow,Mcol,freq)
%        freqs - number vector of spatial frequencies
%        -options-
%        roi - (optional)
%               1 to choose, [w h] to keep size, [left top w h] pos and size
%               Default is the entire image
%        elements - (optional) Which Mueller elements to trace?
%               1 is diagonal. 0 is all. 0 is Default
%        errBars  - (optional) Plot shaded error bars?
%               1 is is to plot error bars. 0 is not. 0 is Default.
%        rawm11   - (optional) Plot m11 not normalized on Trace?
%               1 is to plot it. 0 is not. 0 is Default.
%        f - (optional) figure handle
%
% OUTPUT:
%        roi - roi [left top w h]
%        MNmeans - averaged images -> (4,4,freqs) Normed to M11
%        MNstds  - standard deveiations -> (4,4,freqs) Normed to M11
%        f - figure handle

function varargout = plotMMIvsSI(M,freqs,varargin)
%handle figure handle
if nargin <7, f = figure;
else,      f = varargin{5};
end
figure(f);
if nargout == 4, varargout{4} = f; end
%% Handle ROI
roi = [];
if nargin > 2
    if length(varargin{1}) == 4     %ROI provided
        roi = varargin{1};
    elseif length(varargin{1}) == 2 %ROI dimensions provided
        roiW = varargin{1}(1);
        roiH = varargin{1}(2);
        F = figure; imagesc(M(:,:,1,1,1)); axis equal;
        hr = imrect(gca,[size(M,1)/3 size(M,2)/3  roiW roiH]);
        setResizable(hr,false);
        hm = msgbox({['Drag the rectangle to choose the center of a '...
            num2str(roiW) 'x' num2str(roiH) ' pixel ROI, then'];...
            ['                                    ' ...
            'DOUBLE CLICK it to finish.']},'Help Me Help You');
        uiwait(hm);
        roi = wait(hr); close(F);
    elseif varargin{1} == 1         %ROI must be square
        F = figure; imagesc(M(:,:,1,1,1)); axis equal;
        hr = imrect(gca,[size(M,1)/3 size(M,2)/3 size(M,1)/10 size(M,2)/10]);
        setFixedAspectRatioMode(hr,true);
        msgbox({'Move and resize the rectangle to choose your ROI, then';...
            ['                                    ' ...
            'DOUBLE CLICK it to finish.']},'Help Me Help You');
        roi = wait(hr); close(F);
    end
    %Crop it
    Mc = imcrop(M(:,:,1,1,1),roi);
    Mc = zeros(size(Mc,1), size(Mc,2),4,4,size(M,5));
    for fr = 1:length(freqs)
        for row = 1:4
            for col = 1:4
                Mc(:,:,row,col,fr) = imcrop(M(:,:,row,col,fr),roi);
            end
        end
    end
    
    if nargout >0, varargout{1} = roi; end
else
    Mc = M;
end
%% Plot spatial frequency dependencies

%list of styles for plotting
colors = [0 0 0; .85 .33 .1; .93 .69 .13; .49 .18 .56;...
    .47 .67 .19; .3 .75 .93; .64 .08 .18; 0 1 0];

lstyle = {'-',':'};

%m11 for each spatial frequency
MNmeans = squeeze(mean(mean(Mc,1),2));
m11 = squeeze(MNmeans(1,1,:));

MNstds = zeros(4,4,length(freqs));
for col = 1:4
    for row = 1:4
        for fr = 1:length(freqs)
            if col==1 && row==1
                MNmeans(row,col,fr) = m11(fr);
                MNstds(row,col,fr) = std2(squeeze(Mc(:,:,1,1,fr)));
            else
                MNmeans(row,col,fr) = MNmeans(row,col,fr)./m11(fr);
                MNstds(row,col,fr)  = std2(squeeze(Mc(:,:,row,col,fr))./m11(fr));
            end
        end
    end
end
%standard deviation (of the MEAN)... divide by sqrt(N)
MNstds = MNstds / sqrt(numel(Mc(:,:,1,1,1)));
if nargout>1, varargout{2} = MNmeans; end
if nargout>2, varargout{3} = MNstds; end

%% Check params
if nargin>3
    if varargin{2} %Plot diagonal elements
        PlotDiag = 1;
    else
        PlotDiag = 0;
    end
else
    PlotDiag = 0;
end
if nargin>4
    if varargin{3} == 1  %with stds
        PlotStds = 1;
    else
        PlotStds = 0;
    end
else
    PlotStds = 0;
end
if nargin > 5           %show un-normalized M11?
    if varargin{4}
        plotM11 = 1;
%         m11NormFactor = 2; %this is arbitrary.. for plotting within [-1 1]
    else
        plotM11 = 0;
    end
else
    plotM11 = 0;
end

%% Build Plot
if PlotDiag %Plot diagonal elements
    lineSt = {'-r','-g','-b'};
    y = [squeeze(MNmeans(2,2,:)) squeeze(MNmeans(3,3,:)) squeeze(MNmeans(4,4,:))];
    err = [squeeze(MNstds(2,2,:)) squeeze(MNstds(3,3,:)) squeeze(MNstds(4,4,:))];
    if PlotStds
        if plotM11
            P = shadedErrorBar(freqs,m11./m11(1), MNstds(1,1,:), 'lineprops','k:');
%             P = shadedErrorBar(freqs,m11/m11NormFactor, MNstds(1,1,:), 'lineprops','k:');
            hold on; plotM11 = 1; ylim([-1 1]);
        end
        for L = 1:3
            Ph(L) = shadedErrorBar(freqs, y(:,L), err(:,L), 'lineprops', lineSt{L});
            hold on;
        end
        if plotM11
            legend([P.mainLine,Ph(1).mainLine,Ph(2).mainLine,Ph(3).mainLine],'m11','m22','m33','m44');
        else
            legend([Ph(1).mainLine,Ph(2).mainLine,Ph(3).mainLine],'m22','m33','m44');
        end
    else
        if plotM11
            plot(freqs,m11./m11(1),'k:'); legend('m11'); 
            ylim([-1 1]); hold on;
        end
        plot(freqs,y);
        legend({'m22','m33','m44'});
    end
    legend 'boxoff';
else %Plotting all Elements, no stds
    if PlotStds
        if ~plotM11, MNmeans(1,1,:) = 1; 
        else, MNmeans(1,1,:) = MNmeans(1,1,:)./MNmeans(1,1,1); 
        end
        for col = 1:4
            for row = 1:4
                Ph(row+4*(col-1)) = shadedErrorBar(freqs,squeeze(MNmeans(row,col,:)),...
                    squeeze(MNstds(row,col,:)),'lineprops',...
                    {'LineStyle',lstyle{floor(col/3)+1},...
                    'Color',colors(mod(row+4*(col-1)+7,8)+1,:)}); hold on;
            end
        end
        MNmeans(1,1,:) = m11;
        Pdata = [];
        for i = 1:16, Pdata =[Pdata Ph(i).mainLine]; end
        legend(Pdata,'m11','m21','m31','m41','m12','m22','m32','m42',...
            'm13','m23','m33','m43','m14','m24','m34','m44');
        legend 'boxoff';
    else
        for col = 1:4
            for row = 1:4
                plot(freqs,squeeze(MNmeans(row,col,:)),...
                    'LineStyle',lstyle{floor(col/3)+1},...
                    'Color',colors(mod(row+4*(col-1)+7,8)+1,:)); hold on;
            end
        end
        columnlegend(4,{'m11','m21','m31','m41','m12','m22','m32','m42',...
            'm13','m23','m33','m43','m14','m24','m34','m44'},'boxoff');
    end
end
%or
%     plot(Freqs,reshape(squeeze(mean(mean(M,1),2)),[16 numFreq])./repmat(m11s',[16,1]))
xlabel('spatial frequency (mm^-^1)');
ylabel('element value (au)');


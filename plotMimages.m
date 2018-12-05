%function f = plotMimages(MMI,(NormFlag),(optionPairs))
%
% This function plots the Mueller Matrix of a image data set (4D).
%
% INPUT:
%         MMI(iheight,iwidth,Mrow,Mcol)
%           MMI is a stack with image height 'iheight' and width 'iwidth'.
%           The Mueller matrix row Mrow and column Mcol are assumed to be
%           4x4.
%
%         NormFlag (flag) (optional)
%           If flag 0 is given, no normalization occurs. Raw M  plotted.
%           If flag 1 is given, it will divide all by (1,1) image (besides
%           m11)
%           If flag >1 is given, it does same as above, but divide m11 by
%           this normalization factor.
%           IF flag == 1001, this will norm all elements, and normalize m11
%           by the 95th percentile (to avoid hotspot values)
%
% Option Pairs
%         crange (vector)
%           color range to use for all images: [clow chigh]
%           Defaults to greatest auto-range needed
%
%         cmap (name), standard MATLAB options
%           Defaults to gray (standard montage default)
%
%         fhandle (figure handle)
%           Use given figure handle to put subplots into
%           Defaults to creating a new figure
%
% OUTPUT: f
%           f is the figure handle that holds the plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = plotMimages(MMI,varargin)

%% Parse input and set defaults
p = inputParser;

defaultNormFlag = 0;
numchk    = {'numeric'};
NormAtts  = {'>=',0,'<=',100000};
funName   = 'Rankings';

defaultCrange = [];
CrangeAtts    = {'increasing','numel',2};

defaultCmap = 'parula';
validCmaps  = {'parula','hsv','gray','bone','coper','pink','white','flag',...
    'lines','colorcube','vga','jet','prism','cool','autumn','hot'...
    'spring','winter','summer'};
checkCmap = @(x) any(validatestring(x,validCmaps));

defaultFhandle = [];
checkFhandle = @(x) isa(x,'matlab.ui.Figure');

%validate required input
MMIattributes = {'size',[NaN,NaN,4,4],'ndims', 4};
checkMMI = @(x) validateattributes(x,{'numeric'},MMIattributes);

%ADD scheme
addRequired(p,'MMI',checkMMI);
addOptional(p,'NormFlag',defaultNormFlag,@(x)validateattributes(x,numchk,NormAtts,funName));
addParameter(p,'crange',defaultCrange,@(x)validateattributes(x,numchk,CrangeAtts));
addParameter(p,'cmap',defaultCmap,checkCmap);
addParameter(p,'fhandle',defaultFhandle,checkFhandle);

%Parse input
parse(p,MMI,varargin{:});

%% Build Figure
if isempty(p.Results.fhandle), f = figure;
else,                          f = p.Results.fhandle;
end
figure(f);

if nargout>0, varargout{1} = f; end

%preprocess if flagged for normalization
if p.Results.NormFlag
    %flip factor for flipping m33 and m44, if value is 1111
    if p.Results.NormFlag==1111, flip = ones(4,4); flip(3,3) = -1; flip(4,4)=-1;
    else
        flip = ones(4,4);
    end
    for row = 4:-1:1
        for col = 4:-1:1
            if ~(row==1 && col==1)
                MMI(:,:,row,col) = flip(row,col)*MMI(:,:,row,col)./MMI(:,:,1,1);
            else %special m11 normalizations
                if p.Results.NormFlag == 1001 || p.Results.NormFlag == 1111 %95th percentile norm
                    Mnow = MMI(:,:,1,1);
                    Mnow = sort(Mnow(:),'descend');
                    MMI(:,:,row,col) = MMI(:,:,row,col)./Mnow(round(length(Mnow)*.05+1));
                    
                elseif p.Results.NormFlag == 1002 %95th percentile norm scale [-1 1]
                    Mnow = MMI(:,:,1,1);
                    Mnow = sort(Mnow(:),'descend');
                    MMI(:,:,row,col) = (MMI(:,:,row,col)./Mnow(round(length(Mnow)*.05+1))...
                        -0.5)*2;
                else
                    MMI(:,:,row,col) = MMI(:,:,row,col)./p.Results.NormFlag;
                end
            end
        end
    end
end


%Use Crange or set it
if isempty(p.Results.crange)
    minVal = min(min(min(min(MMI))));
    maxVal = max(max(max(max(MMI(MMI~=1)))));
    if     minVal<0 && abs(minVal)>abs(maxVal), maxVal = -1*minVal;
    elseif minVal<0 && abs(minVal)<abs(maxVal), minVal = -1*maxVal;
    end
    if maxVal > 1, minVal = -1; maxVal = 1; end
    crange = [minVal maxVal];
else
    crange = p.Results.crange;
end


%Build Image Array
%permute matrix so that all vertical components are first
MMP = permute(MMI, [1 3 2 4]);
imagesc(reshape(MMP,[size(MMP,1)*4 size(MMP,3)*4]),crange);
colormap(p.Results.cmap); axis equal; axis off; colorbar;



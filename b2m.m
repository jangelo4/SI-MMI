% function M = b2m(B,A,W)
%
% This function is used to calculate the Mueller matrix M using the
% measurement matrix B and the calibrated PSA and PSG matrices A and W. It
% uses simple matrix multiplication and was created to simplify code.
%
% NOTE: Helpful to use tif2Bmat to create B. Use get_Calib_data then
% System_Calibration_script to get A and W.
%
% INPUTS:
%           B(Iheight,Iwidth,Mrow,Mcol,(f)) This is the measurement matrix 
%           in the format described in tif2Bmat. Effectively a matrix of
%           images.
%
%           A and W (Mrow,Mcol) Point-based calibration given from
%           System_Calibration_script.
%
% OUTPUT:   M(Iheight,Iwidth,Mrow,Mcol) This is the Mueller matrix (of
%           images) for the sample in B.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = b2m(B,A,W)

M = zeros(size(B));
for f = 1:size(B,5)
    for x = 1:size(B,1)
        for y = 1:size(B,2)
            M(x,y,:,:,f) = A\squeeze(B(x,y,:,:,f))/W;
        end
    end
end
function [ stack ] = ReadTIFFobj( PathName, firstframe, lastframe )
%% This function reads TIFF xyz-stacks into a 3D matrix.
% ----------------------------------------------------------------------
% Author: Dr. Knut Kirmse, 2019-2022.
% ----------------------------------------------------------------------
% Input:
% PathName: path & file name of image stack (*.tif)
% firstframe: first frame to be read (optional)
% lastframe: last frame to be read (optional)
% ----------------------------------------------------------------------
% Output:
% stack: image stack
% ----------------------------------------------------------------------
info = imfinfo(PathName);
if nargin < 2
    firstframe = 1;
end
if nargin < 3
    lastframe = numel(info);
end
% ImageJ introduces custom tags causing tifflib to produce warnings.
w = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off', 'MATLAB:imagesci:tifftagsread:expectedTagDataFormat');
obj = Tiff(PathName, 'r');
Nframes = lastframe - firstframe + 1;
stack = zeros(info(1).Height, info(1).Width, Nframes, 'like', obj.read());
count = 1;
for k = firstframe:lastframe
    obj.setDirectory(k);
    stack(:, :, count) = obj.read();
    count = count + 1;
end
obj.close()
warning(w);
end
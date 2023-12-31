function imgsave32(IM,filename,app,ums)
% Save single/double image to a 32bit TIFF file with ImageJ tags
% using a direct interface to libtiff.
% Note that image is appended if filename exists and append condition is 1
%
%   imgsave32(IM, filename)
%
% Input arguments:
%
%   IM         ...   [m x n]   image
%   filename   ...   [string]  file name
%   app        ...    boolean  append condition
%   ums        ...    number   scale in um for ImageJ tag format
%
% See also imgsave16

% Modified from SIMToolbox code:
% Copyright � 2009-2015 Pavel Krizek
%
% This file is part of SIMToolbox.
%
% SIMToolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% SIMToolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with SIMToolbox.  If not, see <http://www.gnu.org/licenses/>.

if isfile(filename) && app==1
    t = Tiff(filename,'a');
else
    t = Tiff(filename,'w');
end

XYRes = 1/ums;

% header
tagstruct.ImageLength     = size(IM,1);
tagstruct.ImageWidth      = size(IM,2);
tagstruct.ResolutionUnit  = Tiff.ResolutionUnit.None;
tagstruct.XResolution     = XYRes;
tagstruct.YResolution     = XYRes;
tagstruct.BitsPerSample   = 32;
tagstruct.SampleFormat    = Tiff.SampleFormat.IEEEFP;
tagstruct.SamplesPerPixel = 1;
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression     = Tiff.Compression.None;
tagstruct.ImageDescription = ['ImageJ=1.53a' newline 'unit=um' newline];

t.setTag(tagstruct)
% write the data
t.write(single(IM));
t.close();

%eof
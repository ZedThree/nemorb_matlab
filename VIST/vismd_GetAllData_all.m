%
% Gets multi-dim data set in HDF5 "CRPP standard format result file".
%
%
% INPUTS: 
%   FileName    : HDF5 file name
%   dim         : = '1d' (resp. '2d') for 1-dim (resp. 2-dim) variables.
%   DataSetName : Name of data set (without path over group hierarchy)
%   group       : = name of species
%
% OUTPUTS:
%   data : Actual data values for all times
%   text : Title string attribute of dataset 
%   PlotOrder : Integer attribute defining plotting order
% 
function [data_full, time, text, PlotOrder,coord1,coord2] = vismd_GetAllData_all(FileName, dim, group, DataSetName);

GroupName = ['/data/var', dim, '/', group, '/', DataSetName];
if nargout > 2,
  % Get title string attribute
  AttrName = [GroupName, '/title'];		 
  attr = hdf5read(FileName, AttrName);
  text = attr.Data;

end

if nargout > 3,
  % Get PlotOrder integer attribute
  AttrName = [GroupName, '/PlotOrder'];		 
  PlotOrder = hdf5read(FileName, AttrName);
end

if nargout > 4,
  AttrName = [GroupName, '/coord1'];		 
  coord1 = hdf5read(FileName, AttrName);
end

if nargout > 5,
  AttrName = [GroupName, '/coord2'];		 
  coord2 = hdf5read(FileName, AttrName);
end
% Read the time
  try
     time = hdf5read(FileName, '/data/var', dim, '/time');
  catch
     s  = lasterror;
     % If there isn\'t a time variable for all the 2d quantities,
     % there might be one in each directory.
     if findstr(s.message,'not an attribute or a dataset')==[]  
       rethrow(s);
     else 
        TimeName = [GroupName, '/time'];
        time = hdf5read(FileName, TimeName);
     end
  end
%Try to read the first frame to determine the size 
try
TimeFrameName = [GroupName, num2str(1, '/%06i')];
data_t = hdf5read(FileName, TimeFrameName);
catch
%Must be single-array format
data_full = double(hdf5read(FileName, [GroupName, '/data']));
return
end
%Read the version number to check whether to do transpose.
k = version;
if ~strcmp(k(1:3),'7.1'),
    data = data_t;
else
    data = transpose(data_t);    
end
if dim=='1d'
data_full = zeros([size(data,1),size(time,2)]);
else
data_full = zeros([size(data),size(time,2)]);
end
% Loop over all the data and collect it into an array.
for i_frame=1:(size(time,2))
TimeFrameName = [GroupName, num2str(i_frame, '/%06i')];

% Earlier versions of matlab transpose the array when they
% read it from the hdf5 file.
if ~strcmp(k(1:3),'7.1'),
    %u = '1 dont transp'
    data = hdf5read(FileName, TimeFrameName);
else
    %u = '1 transp'
    data = transpose(hdf5read(FileName, TimeFrameName));    
end
keyboard

if dim=='2d'
data_full(:,:,i_frame) = data;
else
data_full(:,i_frame) = data;
disp('woop')
end 
end



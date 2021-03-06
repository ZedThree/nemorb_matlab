%
% Gets 0-dim data set in HDF5 "CRPP standard format result file".
%
%  [data, text, PlotOrder] = vis0d_GetDataSet(FileName, DataSetName);
%
% INPUTS: 
%   FileName    : HDF5 file name
%   DataSetName : Name of 0-dim data set (without path over group hierarchy)
%
% OUTPUTS:
%   data : Actual data values for all times
%   text : Title string attribute of dataset 
%   PlotOrder : Integer attribute defining plotting order
% 

function [data, text, PlotOrder] = vis0dgeneric_GetDataSet(FileName, DataSetName)

DataSetName = ['/data/var0d/generic/', DataSetName];

% Get data
data = hdf5read(FileName, DataSetName);

if nargout > 1,
  % Get title string attribute
  AttrName = [DataSetName, '/title'];		 
  attr = hdf5read(FileName, AttrName);
  text = attr.Data;
end

if nargout > 2,
  % Get PlotOrder integer attribute
  AttrName = [DataSetName, '/PlotOrder'];		 
  PlotOrder = hdf5read(FileName, AttrName);
end

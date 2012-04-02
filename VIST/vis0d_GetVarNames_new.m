% Given the group hierarchy of the HDF5 "CRPP standard format result
% file", provides the list of all 0-dim variables
%
% Invoke with 
%	 > list = vis0d_GetVarNames(GroupHierarchy)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%

function list = vis0d_GetVarNames(GroupHierarchy)

% Get list of all data sets in group /data/var0d
list.generic = hdf5_list_all_datasets(GroupHierarchy, '/data/var0d/generic/', 'short');
list.main_ions = hdf5_list_all_datasets(GroupHierarchy, '/data/var0d/main_ions/', 'short');
% Remove time data set from list ;
id_time = strmatch('time', list.generic, 'exact');
list.generic = [list.generic(1:id_time-1), list.generic(id_time+1:end)];

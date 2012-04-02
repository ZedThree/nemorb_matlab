% Given the group hierarchy of the HDF5 "CRPP standard format result
% file", provides the list of all 0-dim variables
%
% Invoke with 
%	 > list = vis0d_GetVarNames(GroupHierarchy, speciesname)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%

function list = vis0d_GetVarNames_all(GroupHierarchy, speciesname)

for kk=1:length(speciesname)
% Pathname for this species
spath=strcat('data/var0d/',speciesname{kk});
% Get list of all data sets in group /data/var0d
list.(speciesname{kk})= hdf5_list_all_datasets(GroupHierarchy, spath, 'short');
end
% Remove time data set from list (assume it's in /data/var0d/generic)
id_time = strmatch('time', list.generic, 'exact');
list.generic = [list.generic(1:id_time-1), list.generic(id_time+1:end)];
% Given the group hierarchy of the HDF5 "CRPP standard format result
% file", provides the list of all multi-dimensional variables
%
% Invoke with 
%	 > list = vismd_GetVarNames(GroupHierarchy, dim)
%
% where the argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%
% and the argument dim is set to '1d' (resp. '2d') if one intends to inquire
% 1-dim (resp. 2-dim) variables.
%

function list = vismd_GetVarNames(GroupHierarchy, dim)

% 1-dim & 2-dim variables are saved their respective groups ...

switch dim,
case '1d',
 id_grp_path = hdf5_identify_group_path(GroupHierarchy, '/data/var1d/');
case '2d'
 id_grp_path = hdf5_identify_group_path(GroupHierarchy, '/data/var2d/');
end

GroupHierarchy = GroupHierarchy.Groups(id_grp_path(1)).Groups(id_grp_path(2));
list.generic = hdf5_list_all_groups(GroupHierarchy.Groups(1));
list.main_ions = hdf5_list_all_groups(GroupHierarchy.Groups(2));

% Remove group '/data/var[1,2]d' from list
list.generic = list.generic(2:end);
list.main_ions = list.main_ions(2:end);				    

% Remove path from names
for iv = 1:length(list.generic),
  is = strfind(list.generic{iv}, '/');
  list.generic{iv} = list.generic{iv}(is(end)+1:end);
end
for iv = 1:length(list.main_ions),
  is = strfind(list.main_ions{iv}, '/');
  list.main_ions{iv} = list.main_ions{iv}(is(end)+1:end);
end

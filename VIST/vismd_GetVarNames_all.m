% Given the group hierarchy of the HDF5 "CRPP standard format result
% file", provides the list of all multi-dimensional variables
%
% Invoke with 
%	 > list = vismd_GetVarNames(GroupHierarchy, dim, speciesname)
%
% where the argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename),
%
% the argument dim is set to '1d' (resp. '2d') if one intends to inquire
% 1-dim (resp. 2-dim) variables, and speciesname is the list of
% kinetic species.
% 
% Output is the list and number of variables stored for each species

function [list nmd] = vismd_GetVarNames_all(GroupHierarchy, dim, speciesname)

% 1-dim & 2-dim variables are saved their respective groups ...

switch dim,
case '1d',
 id_grp_path = hdf5_identify_group_path(GroupHierarchy, '/data/var1d/');
case '2d'
 id_grp_path = hdf5_identify_group_path(GroupHierarchy, '/data/var2d/');
end

GroupHierarchy = GroupHierarchy.Groups(id_grp_path(1)).Groups(id_grp_path(2));

for kk=1:length(speciesname)
list.(speciesname{kk}) = hdf5_list_all_groups(GroupHierarchy.Groups(kk));
% Remove group '/data/var[1,2]d' from list
list.(speciesname{kk}) = list.(speciesname{kk})(2:end);
% Remove path from names
for iv = 1:length(list.(speciesname{kk})),
  is = strfind(list.(speciesname{kk}){iv}, '/');
  list.(speciesname{kk}){iv} = list.(speciesname{kk}){iv}(is(end)+1:end);
end
nmd.(speciesname{kk})=length(list.(speciesname{kk}));
end
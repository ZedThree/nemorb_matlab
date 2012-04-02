%
% Identifies indice id of group name at given level of group hierarchy. 
%
% Invoke with 
% > id = hdf5_identify_group_name(GroupHierarchy, group_name)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%

function id = hdf5_identify_group_name(GroupHierarchy, group_name)

id = -1;
for ig = 1:length(GroupHierarchy),
  if strcmp(GroupHierarchy(ig).Name, group_name),
%    fprintf('Found group name %s at index %i\n', group_name, ig)
    id = ig;
    return
  end
end    
fprintf('Could not find group name %s\n', group_name)

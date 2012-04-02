%
% Identifies array of indices id of group names for all levels of a
% group path.
%
% Invoke with 
% > id = hdf5_identify_group_path(GroupHierarchy, group_path)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%
%

function id = hdf5_identify_group_path(GroupHierarchy, group_path)

group_path = strtrim(group_path);

if group_path(1)   ~= '/', group_path = ['/', group_path] ; end
if group_path(end) ~= '/', group_path(end+1) = '/'; end

is = strfind(group_path, '/');

for ig = 1:length(is)-1,
  group_name = group_path(1:is(ig+1)-1);
  id(ig) = hdf5_identify_group_name(GroupHierarchy.Groups, group_name);
  GroupHierarchy = GroupHierarchy.Groups(id(ig));	
end

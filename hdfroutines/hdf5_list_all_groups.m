%
% Lists the hierarchy of groups of an hdf5 file. 
% This recursive function is invoked as
%
%     list = hdf5_list_all_groups(GroupHierarchy)
%
% The argument GroupHierarchy is the corresponding field of the 
% structure info provided by 
% 
%     info = hdf5info(filename) 
%

function list = hdf5_list_all_groups(GroupHierarchy, list)

%fprintf('%s\n', GroupHierarchy.Name)
if exist('list'),
  list{end+1} = GroupHierarchy.Name;
else
  list{1} = GroupHierarchy.Name;
end

for ig = 1:length(GroupHierarchy.Groups)
list = hdf5_list_all_groups(GroupHierarchy.Groups(ig), list);
end

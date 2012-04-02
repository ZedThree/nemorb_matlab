%
% Given the group path and GroupHierarchy structure of an 
% hdf5 file, povides list of all dataset names in the group.
%
% Invoke with 
% > list = hdf5_list_all_datasets(GroupHierarchy, group_path, style)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%
% OPTIONAL argument style = 'full' -> give dataset names with full
%                                     group paths (default).
%
%                         = 'short' -> drop group path
%

function list = hdf5_list_all_datasets(GroupHierarchy, group_path, style)

if nargin == 2,
  style = 'full';
elseif nargin == 3,
  if ~strcmpi(style, 'full') && ~strcmpi(style, 'short'),
    fprintf('style "%s" not recognized; style must either be "full" or "short"\n', style)
    return
  end
elseif (nargin < 2) && (nargin > 3),
  disp('Not the right number of input arguments')
  return
end

id_grp_path = hdf5_identify_group_path(GroupHierarchy, group_path);

for ig = 1:length(id_grp_path),
  GroupHierarchy = GroupHierarchy.Groups(id_grp_path(ig));
end

for ids = 1:length(GroupHierarchy.Datasets),
  ds_name = GroupHierarchy.Datasets(ids).Name;
  if strcmpi(style, 'full'),
    list{ids} = ds_name;
  elseif strcmpi(style, 'short')
    is = strfind(ds_name, '/');
    list{ids} = ds_name(is(end)+1:end);
  end
end

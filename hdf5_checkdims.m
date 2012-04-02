function dims = hdf5_checkdims(info, datasetname, ind)

% Matlab has a tendancy to crash if it tries to read a NULL hdf5
% array or dataset, so this routine checks the size of
% datasetname. If dims=0, then array is probably NULL
%
% Inputs:
%   info        : HDF5 file info
%   datasetname : Dataset to check - give full path
%                 (e.g. '/data/var2d/generic/efspecmn/time')
%   ind         : array containing indices of simulations, defaults
%                 to 1
% Output:
%   dims        : size of array
%
% To check non-scalar variables, check their '/time' dataset

if (~exist('ind'))
ind=1;
end

for i=1:length(ind) % loop over ind
k=ind(i);

% find indices of slashes in datasetname
temp=findstr('/',datasetname);
% add on index of datasetname(end)
temp(end+1)=length(datasetname)+1;

GroupHierarchy=info(k).GroupHierarchy;

% loop over number of groups in datasetname
for ii=1:length(temp)-2
% names of groups are located between adjacent slashes
groupname=datasetname(1:temp(ii+1)-1);

% get index of current group name
for jj=1:length(GroupHierarchy.Groups)
if strcmp(GroupHierarchy.Groups(jj).Name, groupname) == 1
  igroup=jj;
end
end
% move down into groupname
GroupHierarchy=GroupHierarchy.Groups(igroup);
end % group loop

for jj=1:length(GroupHierarchy.Datasets)
if strcmp(GroupHierarchy.Datasets(jj).Name, datasetname) ==1
  iset=jj;
end
end

% check number of dimensions
dims=GroupHierarchy.Datasets(iset).Dims;

end % loop on ind
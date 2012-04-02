% Given the group hierarchy of the HDF5 "CRPP standard format result
% file", provides the list of all species names in a given directory
%
% Invoke with 
%	 > speciesname = GetSpeciesNames(GroupHierarchy, dir)
%
% The argument GroupHierarchy is the corresponding field of the
% structure info provided by
%   
%      info = hdf5info(filename)
%

function speciesname = GetSpeciesNames(GroupHierarchy, dir)

for ii=1:length(GroupHierarchy.Groups)
%see if the name of the group i matches "data"
if strcmp(GroupHierarchy.Groups(ii).Name, dir) == 1
   idir=ii;
end
end

% length of data dir name
dirnamel=length(GroupHierarchy.Groups(idir).Groups(1).Name);

% Loop over species
for kk=1:length(GroupHierarchy.Groups(idir).Groups(1).Groups)
% Get name of species and strip off the directory name
speciesname(kk)={GroupHierarchy.Groups(idir).Groups(1).Groups(kk).Name}
speciesname{kk}=speciesname{kk}(dirnamel+2:end)
end
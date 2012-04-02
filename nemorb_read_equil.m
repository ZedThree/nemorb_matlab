function [nemo_t]=nemorb_read_equil(nemo_t, info, ind)
%######################################
%[nemo_t]=nemorb_read_equil(nemo_t, info, ind)
%######################################
%---------------
%
%Description
%
% Read the various equil quantities 
%
%---------------
%
%Input arguments
%
% nemo_t = array containing simulation data structure
% info   = array containing simulation info structure
% ind    = array containing indices of simulations to read equil data from
%---------------
%
%Output
%
% nemo_t = simulation data structure with equil data
%
%---------------

global timelab
if exist('ind') == 0
  ind=1;
end

pwd_old=pwd;
for i=1:length(ind) % Loop on simulations
  k=ind(i);
  cd(nemo_t(k).path);

  disp('reading equilibrium profiles')
  %find the index of the group "equil"
  for ii=1:length(info(k).GroupHierarchy.Groups)
    %see if the name of the group i matches "equil"
    if strcmp(info(k).GroupHierarchy.Groups(ii).Name,'/equil') == 1
      iequil=ii;
    end
  end

  %find the index of the group "equil/sc"
  for ii=1:length(info(k).GroupHierarchy.Groups(iequil).Groups)
    %see if the name of the group i matches "/equil/sc"
    if strcmp(info(k).GroupHierarchy.Groups(iequil).Groups(ii).Name,'/equil/sc')==1
      isc=ii;
    end
  end

  % Length of profiles dir name
  dirnamel=length(info(k).GroupHierarchy.Groups(iequil).Groups(isc).Name);

  %Number of profiles for this species
  nprof=size(info(k).GroupHierarchy.Groups(iequil).Groups(isc).Datasets,2);

  for ii=1:nprof
    %location of equil/profiles in hdf5 file
    pathname=info(k).GroupHierarchy.Groups(iequil).Groups(isc).Datasets(ii).Name;
    %Check dataset isn't empty - matlab doesn't like it
    dims=hdf5_checkdims(info,pathname,k);
    if dims==0
      err=sprintf('%s','ERROR: ', pathname, ' is empty, skipping');
      disp(err)
      continue
    end
    %read value of input parameter
    varvalue=hdf5read(nemo_t(k).filename,pathname);
    pathl=length(pathname);
    %remove '/equil/profiles/species-name/' from the pathname
    listprof{ii}=pathname(dirnamel+2:pathl);
    %check if variable already exists
    if exist(listprof{ii})==1
      sprintf('%s','ERROR: ', speciesname(kk).(listprof{ii}),' already exists')
      i
      return
    end
    % Create structure attributes
    % Brackets around speciesname{kk} and listprof{ii} very important!
    % Let's us use the string in the structure name
    nemo_t(k).equil.(listprof{ii})=varvalue;
    clear varvalue
  end


end % Loop on sims
nemo_t=nemo_t;
cd(pwd_old);

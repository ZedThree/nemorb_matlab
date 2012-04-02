function [struct_out, info]=nemorb_load(path, name, mydim)
%######################################
%[struct_out, info]=nemorb_load(path, name, mydim)
%######################################
%---------------
%
%Description
%
%Routines that load hdf5 result files of simulations named through array name and defined through directories named path 
%
%---------------
%
%Dummy arguments
%
%path = array containing the directories where the hdf5 files are stored
%name = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%struct_out = structure containing 0D data and 1D data information
%info       = structure containing the hierarchy of the hdf5 file
%---------------
%EG: 
%    [nemo,info]=nemorb({'.'},{'nemo'})       reads 0d and 1d data
%    [nemo,info]=nemorb({'.'},{'nemo'},'0d')  reads 0d data only
%    [nemo,info]=nemorb({'.'},{'nemo'},'1d')  reads 1d data only
%    [nemo,info]=nemorb({'.'},{'nemo'},'no')  reads no data 
%    [nemo,info]=nemorb({'.','./NEMORB/WK/quick/'},{'nemo','quick'})


global timelab

load_0d=1;
load_1d=1;

if exist('mydim')==0
  mydim='undef';
else
  if (mydim=='0d')
    load_1d=0;
  end
  if (mydim=='1d')
    load_0d=0;
  end
  if (mydim=='no')
    load_0d=0;
    load_1d=0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD HDF5 DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn off deprecated syntax warning for hdf5read
warning('off', 'MATLAB:hdf5readc:deprecatedAttributeSyntax');
warning('off','MATLAB:imagesci:hdf5readc:deprecatedAttributeSyntax');
warning('off', 'nemorb:emptydata');

%set default values for graphical output
set(0,'DefaultAxesFontsize',16)
set(0,'DefaultTextFontsize',16)
set(0,'DefaultAxesFontname','helvetica')
set(0,'DefaultLineLinewidth',2.0)
set(0,'DefaultLineColor','black')
set(0,'DefaultLineMarkersize',8)
set(0,'DefaultFigureUnits','centimeters')
set(0,'DefaultFigurePosition',[2 2 24 18])
set(0,'DefaultFigurePaperpositionmode','auto')

set(0,'DefaultAxesFontsize',20)
set(0,'DefaultTextFontsize',20)
set(0,'DefaultAxesFontname','helvetica')
set(0,'DefaultLineLinewidth',3.0)
set(0,'DefaultLineColor','black')
set(0,'DefaultLineMarkersize',8)
set(0,'DefaultFigureUnits','centimeters')
set(0,'DefaultFigurePosition',[2 2 24 18])

%axis coordinates
ax=3.1199;
ay=1.4791;
aw=16.7661;
ah=14.6636;

%Graphical output
timelab{1}='t [\Omega_i^{-1}]';

%name of hdf5 file
filename='orb5_res.h5';
disp(['hdf5 file: ', filename])
pwd_old=pwd;
for jj=1:length(path) %loop on simulation

%Go to simulation directory
cd(path{jj});
disp(['location: ', path{jj}])
%Create first main elements of the structure
struct_out(jj).name=name{jj};
struct_out(jj).path=path{jj};
struct_out(jj).filename=filename;

%Get general information of database
info(jj)=hdf5info(filename);


%%%%%%%%%%%%%%%%%%%DIAGNOSTICS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%starting time for linear growth rate
struct_out(jj).tmin=80000;
%Ending time for linear growth rate
struct_out(jj).tmax=120000;
%starting time for real frequency
struct_out(jj).tmin_rf=5000;
%Ending time for real frequency
struct_out(jj).tmax_rf=80000;
%tension for splines smoothing
struct_out(jj).tol=1e-5;
%dimension of higher resolution grid for 1D profiles
struct_out(jj).nhr=128;
%lower value of s for RH tests
struct_out(jj).smin_RH=0.52;
%upper value of s for RH tests
struct_out(jj).smax_RH=0.72;
%starting time for decay, RH tests
struct_out(jj).tmin_RH=1000;
%ending time for decay, RH tests
struct_out(jj).tmax_RH=14000;
%Normalise chi with Ln (1) or with minor radius a (0)
struct_out(jj).norm_Dimits=0;
%half width of output frequency spectrum
struct_out(jj).wfs=20;
%\chi vs R/LT is computed every tlength time steps
struct_out(jj).tlength=4;
%radial coordinates for plots of transport quantities (1:s, 2: sqrt(V/V_tot)=rho/a for ad hoc)
struct_out(jj).irad_coord=2;
%Minimum rad for profile averaging
struct_out(jj).radlow=0.4;
%Maximum rad for profile averaging
struct_out(jj).radup=0.6;
%lower value of s or rho/a for 1D profiles
struct_out(jj).radmin_1D=0.0;
%upper value of s or rho/a for 1D profiles
struct_out(jj).radmax_1D=1.0;
%lower value of s or rho/a for 1D ZF profiles
struct_out(jj).smin_1D=0.0;
%upper value of s or rho/a for 1D ZF profiles
struct_out(jj).smax_1D=1.0;
%number of isosurfaces for binning of electrons in (s,theta,lambda)
struct_out(jj).nsurf_3D=64;
%number of surfaces for binning of electrons in (theta,lambda)
struct_out(jj).nsurf_2D=4;
%time window for moving average, in a/c_s units
struct_out(jj).tav=400;
%starting time for moving average in omega_i-1 units
struct_out(jj).tav_start=75000;


%%%%%%%%%%%%%%%%%%%%%  GET PARAMETERS OF THE SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp('reading parameters of simulations')
%find the index of the group "parameters"
for ii=1:length(info(jj).GroupHierarchy.Groups)
%see if the name of the group i matches "parameters"
   if strcmp(info(jj).GroupHierarchy.Groups(ii).Name,'/parameters') == 1
      ipara=ii;
      struct_out(jj).ipara=ipara;
   end
end

% length of parameters dir name
paradirnamel=length(info(jj).GroupHierarchy.Groups(ipara).Name);
% Loop over species
for kk=1:length(info(jj).GroupHierarchy.Groups(ipara).Groups)
  % Get name of species and strip off '/parameters/'
  speciesname(kk)={info(jj).GroupHierarchy.Groups(ipara).Groups(kk).Name};
  speciesname{kk}=speciesname{kk}(paradirnamel+2:end);
  %Number of input parameters for this species
  nparam=size(info(jj).GroupHierarchy.Groups(ipara).Groups(kk).Datasets,2);
  %/parameters/species-name/
  dirname=strcat(info(jj).GroupHierarchy.Groups(ipara).Groups(kk).Name,'/');
  dirnamel=length(dirname);

  for ii=1:nparam
    %location of input parameter in hdf5 file
    pathname=info(jj).GroupHierarchy.Groups(ipara).Groups(kk).Datasets(ii).Name;
    %Check dataset isn't empty - matlab doesn't like it
    dims=hdf5_checkdims(info,pathname,jj);
    if dims==0
      warning('nemorb:emptydata','%s is empty, skipping 11', pathname)
      err=sprintf('%s','ERROR: ', pathname, ' is empty, skipping');
      disp(err)
      continue
    end
    %read value of input parameter
    varvalue=hdf5read(filename,pathname);
    pathl=length(pathname);
    %remove '/parameters/species-name/' from the pathname
    listpara{ii}=pathname(dirnamel+1:pathl);
    %check if variable already exists
    if exist(listpara{ii})==1
      sprintf('%s','ERROR: ', speciesname(kk).(listpara{ii}),' already exists')
      i
      return
    end
    kk; % de-mute for testing
    ii;
    % Create structure attributes
    % Brackets around speciesname{kk} and listpara{ii} very important!
    % Let's us use the string in the structure name
    struct_out(jj).(speciesname{kk}).(listpara{ii})=varvalue(1);
    clear varvalue
  end
  clear listpara
end
clear speciesname

% Older files may have data just in /parameters/
% $$$ npara=size(info(jj).GroupHierarchy.Groups(ipara).Datasets,2);
% $$$ dirname=strcat(info(jj).GroupHierarchy.Groups(ipara).Name,'/');
% $$$ dirnamel=length(dirname);
% $$$ for ii=1:nparam
% $$$ %location of input parameter in hdf5 file
% $$$ pathname=info(jj).GroupHierarchy.Groups(ipara).Datasets(ii).Name;
% $$$ %read value of input parameter
% $$$ varvalue=hdf5read(filename,pathname);
% $$$ pathl=length(pathname);
% $$$ %remove '/parameters' from the pathname
% $$$ listpara{ii}=pathname(dirnamel+1:pathl);
% $$$ %check if variable already exists
% $$$ if exist(listpara{ii})==1
% $$$ sprintf('%s','ERROR: ', listpara{ii},' already exists')
% $$$ return
% $$$ end
% $$$ % Create structure attributes
% $$$ struct_out(jj).parameters.(listpara{ii})=varvalue;
% $$$ clear varvalue
% $$$ end
clear speciesname
%%%%%%%%%%%%%%%%%%%%%%%%% GET ANALYTICAL PROFILES %%%%%%%%%%%%%%%%%%%%%%%%%
%
disp('reading equilibrium profiles')
%find the index of the group "equil"
for ii=1:length(info(jj).GroupHierarchy.Groups)
  %see if the name of the group i matches "equil"
  if strcmp(info(jj).GroupHierarchy.Groups(ii).Name,'/equil') == 1
    iequil=ii;
  end
end

%find the index of the group "equil/profiles"
for ii=1:length(info(jj).GroupHierarchy.Groups(iequil).Groups)
  %see if the name of the group i matches "/equil/profiles"
  if strcmp(info(jj).GroupHierarchy.Groups(iequil).Groups(ii).Name,'/equil/profiles')==1
    iprof=ii;
  end
end
% Length of profiles dir name
profdirnamel=length(info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Name);

% Loop over species
for kk=1:length(info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Groups)
  % Get name of species and strip off '/equil/profiles/'
  speciesname(kk)={info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Groups(kk).Name};
  speciesname{kk}=speciesname{kk}(profdirnamel+2:end);
  struct_out(jj).allspecies{kk}=speciesname{kk};
  %Number of profiles for this species
  nprof=size(info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Groups(kk).Datasets,2);
  %/equil/profiles/speciesname
  dirname=strcat(info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Groups(kk).Name,'/');
  dirnamel=length(dirname);

  for ii=1:nprof
    %location of equil/profiles in hdf5 file
    pathname=info(jj).GroupHierarchy.Groups(iequil).Groups(iprof).Groups(kk).Datasets(ii).Name;
    %Check dataset isn't empty - matlab doesn't like it
    dims=hdf5_checkdims(info,pathname,jj);
    if dims==0
      err=sprintf('%s','ERROR: ', pathname, ' is empty, skipping');
      disp(err)
      continue
    end
    %read value of input parameter
    varvalue=hdf5read(filename,pathname);
    pathl=length(pathname);
    %remove '/equil/profiles/species-name/' from the pathname
    listprof{ii}=pathname(dirnamel+1:pathl);
    %check if variable already exists
    if exist(listprof{ii})==1
      sprintf('%s','ERROR: ', speciesname(kk).(listprof{ii}),' already exists')
      i
      return
    end
    kk; % de-mute for testing
    ii;
    % Create structure attributes
    % Brackets around speciesname{kk} and listprof{ii} very important!
    % Let's us use the string in the structure name
    struct_out(jj).(speciesname{kk}).(listprof{ii})=varvalue;
    clear varvalue
  end
  clear listprof
end
clear speciesname

%%%%%%%%%%%%%%%%%%%%%%%%%%% GET EQUILIBRIUM SCALARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%find the index of the group "equil/scalars"
for ii=1:length(info(jj).GroupHierarchy.Groups(iequil).Groups)
  %see if the name of the group i matches "/equil/scalars"
  if strcmp(info(jj).GroupHierarchy.Groups(iequil).Groups(ii).Name,'/equil/scalars') == 1
    iscal=ii;
  end
end

% Length of profiles dir name
scaldirnamel=length(info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Name);

disp('reading equilibrium scalars')
for kk=1:length(info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Groups)
  % Get name of species and strip off '/equil/scalars/'
  speciesname(kk)={info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Groups(kk).Name};
  speciesname{kk}=speciesname{kk}(scaldirnamel+2:end);
  
  %Number of scalars
  nscal=size(info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Groups(kk).Datasets,2);
  %/equil/scalars/species-name/
  dirname=strcat(info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Groups(kk).Name,'/');
  dirnamel=length(dirname);
  for ii=1:nscal
    %location of scalars in hdf5 file
    pathname=info(jj).GroupHierarchy.Groups(iequil).Groups(iscal).Groups(kk).Datasets(ii).Name;
    % Check dataset isn't empty - matlab doesn't like it
    dims=hdf5_checkdims(info,pathname,jj);
    if dims==0
      err=sprintf('%s','ERROR: ', pathname, ' is empty, skipping');
      disp(err)
      continue
    end
    varvalue=hdf5read(filename,pathname);
    pathl=length(pathname);
    %remove '/equil/scalars' from the pathname
    listscal{ii}=pathname(dirnamel+1:pathl);
    %check if variable already exists  !!!!!!!CHECK!!!!!!!!!
    if exist((listscal{ii}))==1
      sprintf('%s','ERROR: ', listscal{ii},' already exists')
      i
      return
    end
    % Create structure attributes
    struct_out(jj).(speciesname{kk}).(listscal{ii})=varvalue(1);
    clear varvalue
  end
  clear listscal
end
clear speciesname

struct_out(jj).generic.aspect_ratio(1)=struct_out(jj).generic.r0_mid(1)/struct_out(jj).generic.a_mid(1)*struct_out(jj).generic.d_norm(1);

if struct_out(jj).generic.aspect_ratio > 100
  disp('WARNING ASPECT RATIO MODIFIED, DUE TO CHANGE OF CODE')
  struct_out(jj).generic.aspect_ratio(1)=struct_out(jj).generic.aspect_ratio(1)/ ...
      struct_out(jj).generic.d_norm(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET 0D DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (load_0d==1)
  
  disp('reading 0D data')
  % May have different species in /data/ than in parameters, so get
  % names again. Look in /data/var0d/ for names
  for ii=1:length(info(jj).GroupHierarchy.Groups)
    %see if the name of the group i matches "data"
    if strcmp(info(jj).GroupHierarchy.Groups(ii).Name,'/data') == 1
      idata=ii;
    end
  end
  % length of data dir name
  datadirnamel=length(info(jj).GroupHierarchy.Groups(idata).Groups(1).Name);
  for kk=1:length(info(jj).GroupHierarchy.Groups(idata).Groups(1).Groups)
    % Get name of species and strip off '/data/vis0d/'
    speciesname(kk)={info(jj).GroupHierarchy.Groups(idata).Groups(1).Groups(kk).Name};
    speciesname{kk}=speciesname{kk}(datadirnamel+2:end);
  end
  
  %Get list of 0d variable name (without the time)
  list0d=vis0d_GetVarNames_all(info(jj).GroupHierarchy, speciesname);
  
  %Get time
  struct_out(jj).time = vis0dgeneric_GetDataSet(filename,'time');
  %Number of 0D diagnostics
  for kk=1:length(speciesname)
    n0d.(speciesname{kk})=length(list0d.(speciesname{kk}));
  end
  
  for kk=1:length(speciesname)
    %Create Matlab variables - generic
    for ii=1:n0d.(speciesname{kk})
      %read data
      %Check dataset isn't empty - matlab will crash if attempts to load
      %an empty dataset
      dims=info(jj).GroupHierarchy.Groups(idata).Groups(1).Groups(kk).Datasets(ii).Dims;
      if dims==0
	warning('nemorb:emptydata','%s is empty, skipping', [speciesname{kk} ...
		    '.' list0d.(speciesname{kk}){ii}])
% $$$   err=sprintf('%s','ERROR: ', speciesname{kk}, '.', ...
% $$$ 	    list0d.(speciesname{kk}){ii}, ' is empty, skipping');
% $$$   disp(err)
	continue
      end
      varvalue = vis0d_GetDataSet_all(filename, ...
				      list0d.(speciesname{kk}){ii},speciesname{kk});
      %check is variable already exists
      if exist(list0d.(speciesname{kk}){ii})==1
	sprintf('%s','ERROR: ', list0d.(speciesname{kk}){ii},' already exists')
	i
	return
      end
      % Create structure attributes
      struct_out(jj).(speciesname{kk}).(list0d.(speciesname{kk}){ii})=single(varvalue);
      clear varvalue
    end %data loop
  end %species loop
  
end %load_0D loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET 1D DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (load_1d==1)
  disp('reading 1D data')
  %1d variables
  [list1d n1d]=vismd_GetVarNames_all(info(jj).GroupHierarchy,'1d', speciesname);
  
  
  for kk=1:length(speciesname)
    for ii=1:n1d.(speciesname{kk})
      % Check dataset isn't empty - matlab will crash if attempts to load
      % an empty dataset. For 1D/2D check "time" dataset (3/4 respectively)
      dims=info(jj).GroupHierarchy.Groups(idata).Groups(2).Groups(kk).Groups(ii).Datasets(3).Dims;
      if dims==0
	warning('nemorb:emptydata','%s is empty, skipping', [speciesname{kk} ...
		    '.' list1d.(speciesname{kk}){ii}])
% $$$ err=sprintf('%s','ERROR: ', speciesname{kk}, '.', ...
% $$$ 	    list1d.(speciesname{kk}){ii}, ' is empty, skipping');
% $$$ disp(err)
	continue
      end
      %read data
      [data_full, time_1D, text, PlotOrder,coord1D]=vismd_GetAllData_all(filename,'1d',speciesname{kk},list1d.(speciesname{kk}){ii});
      %check if variable already exists
      if exist(char(list1d.(speciesname{kk}){ii}))==1
	warning('nemorb:alreadyexists', '%s already exists', list1d.(speciesname{kk}){ii})
	return
      end
      struct_out(jj).(speciesname{kk}).(list1d.(speciesname{kk}){ii})=single(data_full);
      clear data_full
      %assign various 1D grid
      if strcmp('phibar',list1d.(speciesname{kk}){ii}) == 1
	struct_out(jj).sgrid = coord1D;
	struct_out(jj).dsgrid = struct_out(jj).sgrid(2) - ...
	    struct_out(jj).sgrid(1);
      elseif strcmp('efluxw_rad',list1d.(speciesname{kk}){ii})==1
	struct_out(jj).psi = coord1D;
	struct_out(jj).dpsigrid = struct_out(jj).psi(2) - ...
	    struct_out(jj).psi(1);
	% efluxw_rad only has data every (NMAX_SPECIES)-th row
	% This bug has now been fixed, but may still be present in older
	% files.
	% First, check if array is same length as time - if it's some
	% multiple of time, then this multiple tells you how far apart the
	% data is actually spaced.f
% $$$ bugtest=size(struct_out(jj).(speciesname{kk}).efluxw_rad,2)/ ...
% $$$ 	(length(struct_out(jj).time) - 1);
% $$$ if bugtest ~= 1
% $$$ err=sprintf('ERROR: efluxw_rad bug, fixing...');
% $$$ disp(err)
% $$$ idx = 1:bugtest:length(struct_out(jj).(speciesname{kk}).efluxw_rad);
% $$$ struct_out(jj).(speciesname{kk}).efluxw_rad = struct_out(jj).(speciesname{kk}).efluxw_rad(:,idx);
% $$$ end
      elseif strcmp('phi0_chi',list1d.(speciesname{kk}){ii})==1
	struct_out(jj).chi_1D=coord1D;
      elseif strcmp('f_av',list1d.(speciesname{kk}){ii})==1
	struct_out(jj).time_1D=time_1D;
      end % 1D grids
    end % variable loop
  end % species loop
end % load_1d

% Useful to know what species are in data
for ii=1:length(speciesname)
  struct_out(jj).kineticspecies{ii}=speciesname{ii};
end

%
% Get general information on 2D data
%
%2d variables
[list2d n2d]=vismd_GetVarNames_all(info(jj).GroupHierarchy,'2d',speciesname);

struct_out(jj).nsteps=length(struct_out(jj).time);
%struct_out(jj)=orderfields(struct_out(jj));

%if exist('simu')==0
%is=1
%else
%is=length(simu)+1
%end
%eval(['simu(' num2str(is) ')=orderfields(struct_out(jj));']);
%clear all variables  except simu
%clear -regexp ^[^s]* ^s[^i]\w+ ^s[i][^m]\w+ ^s[i][m][^u]
%end %loop on simulations

end %loop on simulations

cd(pwd_old);

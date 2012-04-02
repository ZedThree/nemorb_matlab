function [sim hh] = nemorb_plot_flux(sim, species, flux, ind)
%################################
% [sim hh] = nemorb_plot_flux(sim, species, flux, ind)
%################################
%---------------
%
% Description
%
% Plot flux-sufaced averaged radial flux on (time, s)
%
%---------------
%
% Input arguments
%
% sim     = simulation data structure, created with nemorb_load
% species = name of species to plot flux of
% flux    = name of flux to be plotted
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
% Output
%
% sim     = simulation data structure
% hh      = (optional) figure handle. If present, plots figure in
%           background. Display with set(hh,'Visible','on')
%---------------

global timelab

% Check if we are looping over multiple simulations
if exist('ind')==0
  ind=1;
end

% Remember where we are
pwd_old = pwd;

% Loop over simulations
for ii=1:length(ind)
  k=ind(ii);
  
  % If given output argument hh, plot figure in background
  if exist('hh')==1
    h(k) = figure('Visible','off');
  else
    h(k) = figure;
  end

  % Check flux exists
  if isfield(sim(k).(species),flux)==0
    warning('nemorb:emptydata',['%s does not exist for species %s in simulation' ...
	' %s'], flux, species, sim(k).name)
    err = sprintf('%s does not exist for species %s in simulation %s', flux, ...
	species, sim(k).name)
    disp(err)
    return
  end

  % Need 1D profiles - generate them if necessary
  if isfield(sim(k).(species),'time_1D')==0
    cd(sim(k).path);
    filename = sim(k).filename;
    if sim(k).generic.nsel_equil == 2 & sim(k).(species).nsel_profile > 1
      s1=2;
      s2=1;
    else
      s1=1;
      s2=0;
    end
    sim(k).(species).time_1D = double(hdf5read(filename,['/data/var1d/',species,'/',flux,'/time']));
    tmp = double(hdf5read(filename,['/data/var1d/',species,'/',flux,'/coord1']));
    sim(k).(species).psi_prof1D = tmp(s1:end-s2);
    sim(k).(species).s_prof1D   = sqrt(sim(k).(species).psi_prof1D);
  end
  
  % Plot heat flux
  pcolor(sim(k).(species).time_1D,sim(k).(species).s_prof1D, ...
       double(sim(k).(species).(flux)));
  shading interp;
  colorbar
  ylabel('s')
  xlabel('Time [\Omega_{ci}]')
  title([sim(k).name ' - ' flux])
end

hh=h;
% Make sure 1D profiles get saved
sim=sim;
% Return to whence we came
cd(pwd_old);
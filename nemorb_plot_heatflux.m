function [sim hh] = nemorb_plot_heatflux(sim, species, ind)
%################################
% [sim hh] = nemorb_plot_heatflux(sim, species, ind)
%################################
%---------------
%
% Description
%
% Plot flux-sufaced averaged radial heat flux against time, s
%
%---------------
%
% Input arguments
%
% sim     = simulation data structure, created with nemorb_load
% species = name of species to plot heat flux of
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
  kk=ind(ii);
  
  % If given output argument hh, plot figure in background
  if exist('hh')==1
    hh = figure('Visible','off');
  else
    hh = figure;
  end

  % Need 1D profiles - generate them if necessary
  if isfield(sim(kk).(species),'time_1D')==0
    cd(sim(k).path);
    filename = sim(k).filename;
    if sim(k).generic.nsel_equil == 2 & sim(k).(species).nsel_profile > 1
      s1=2;
      s2=1;
    else
      s1=1;
      s2=0;
    end
    sim(k).(species).time_1D = double(hdf5read(filename,['/data/var1d/',species,'/f_av/time']));
    tmp = double(hdf5read(filename,['/data/var1d/',species,'/efluxw_rad/coord1']));
    sim(k).(species).psi_prof1D = tmp(s1:end-s2);
    sim(k).(species).s_prof1D   = sqrt(sim(k).(species).psi_prof1D);
  end
  
  % Plot heat flux
  pcolor(sim(kk).(species).time_1D,sim(kk).(species).s_prof1D, ...
       double(sim(kk).deuterium.efluxw_rad));
  shading interp;
  colorbar
  ylabel('s')
  xlabel('Time [\Omega_{ci}]')
  title(sim(kk).name)
end

% Make sure 1D profiles get saved
sim=sim;
% Return to whence we came
cd(pwd_old);
function [hh] = nemorb_snr(sim, species, ind)
%################################
% [hh] = nemorb_snr(sim, species, ind)
%################################
%---------------
%
% Description
%
% Plot signal-to-noise ratio (full and non-zonal components) of simulation
%
%---------------
%
% Input arguments
%
% sim     = simulation data structure, created with nemorb_load
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
% Output
%
% hh      = (optional) figure handle. If present, plots figure in
%           background. Display with set(hh,'Visible','on')
%---------------

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
  plotname = [sim(k).name ' - SNR'];
  if exist('hh')==1
    h(k) = figure('Visible','off','Name',plotname);
  else
    h(k) = figure('Name',plotname);
  end

  % Plot SNR
  plot(sim(k).time, sim(k).(species).signal./sim(k).(species).noise)
  hold on
  plot(sim(k).time, sim(k).(species).signal_nonzonal./ ...
       sim(k).(species).noise,'--r')
  plot([sim(k).time(1) sim(k).time(end)],[10 10],':k')
  ylabel('SNR')
  xlabel('Time [\Omega_{ci}]')
  title([sim(k).name ' -  signal-to-noise ratio'])
  legend('Full','Non-zonal','SNR=10','Location','NorthWest')
end

hh=h;
% Return to whence we came
cd(pwd_old);
function sim = nemorb_coarsegrain_diags(sim, species, graphs, ind)
%################################
% sim = nemorb_coarsegrain_diags(sim, species, graphs, ind)
%################################
%---------------
%
% Description
%
% Get various diagnostics of the coarse-graining
%
%---------------
%
% Input arguments
%
% sim     = array containing the simulations data, created with nemorb_load.m
% species = name of species
% graphs  = [ 1 : 0 ] graphs on or off
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
% Output
%
% sim         
%---------------

if exist('ind') == 0
  ind = 1;
end

pwd_old = pwd;

for i=1:length(ind) % Loop on simulations
  k = ind(i);
  cd(sim(k).path);
  
  % Read diagnostic data
  non_cg_markers = hdf5read(sim.filename,['/data/var0d/' species '/cum_non_cg_markers']);
  tot_cg_markers = hdf5read(sim.filename,['/data/var0d/' species '/tot_cg_markers']); 
  tot_cg_bins = hdf5read(sim.filename,['/data/var0d/' species '/tot_cg_bins']);

  % Normalisations - number of markers and bins
  nptot = sim(k).(species).nptot;
  bins_tot = sim.generic.nphi*sim.(species).nbins_en*sim.(species).nbins_lambda;
  bins_tot = bins_tot*sim.generic.nbin_chi_coarse*sim.generic.ns;

  % Find non-zero data points
  % This is because the diagnostic is output at the same frequency as the 0d diagnostics,
  % which is not necessarily the same as the coarse-graining frequency
  ii = find(tot_cg_markers);
  
  if graphs==1
    figure;
    hold on;
    xlabel('Time [\Omega_{ci}^{-1}]')
    plot(sim.time(ii), non_cg_markers(ii)/nptot, 'r')
    plot(sim.time(ii), tot_cg_markers(ii)/nptot, 'k')
    plot(sim.time(ii), tot_cg_bins(ii)/bins_tot, 'g')
    legend('Fraction of non-CG markers so far', 'Fraction of CG markers at time', 'Fraction of CG bins at time')
  end
end

cd(pwd_old);
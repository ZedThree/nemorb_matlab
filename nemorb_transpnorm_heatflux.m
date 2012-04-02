function sim = nemorb_transpnorm_heatflux(sim,transp,species,ind)
%################################
%sim = nemorb_transpnorm_heatflux(sim, transp, species, times, ind)
%################################
%---------------
%
% Description
%
% Normalise heat flux according to TRANSP equilibrium data
%
%---------------
%
% Input arguments
%
% sim     = structure containing the simulations data, created with nemorb_load.m
% transp  = structure containing TRANSP data
% species = name of species to normalise
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
% Output
%
% sim     = the simulation data structure 
%---------------

% SI units
ElectronCharge=1.6022E-19;  % Coulomb

% Check if we are looping over multiple simulations
if exist('ind') == 0
  ind=1;
end

% Loop over simulations
for i=1:length(ind)
  k=ind(i);
  
  % Get actual values of equilibrium quantities at reference surface
  n0   	= spline(transp.s, transp.net, sim(k).(species).speak);
  te_s0	= spline(transp.s, transp.tet, sim(k).(species).speak);
  cs0  	= spline(transp.s, transp.cs,  sim(k).(species).speak);

  % Do normalisation
  sim(k).(species).q_rad=sim(k).(species).efluxw_rad*n0*te_s0*cs0*ElectronCharge;
  % For electrons, also normalise trapped+passing
  if strcmp(species,'electrons')
    efluxw_rad_el=sim(k).(species).efluxw_rad_tr+sim(k).(species).efluxw_rad_p;
    sim(k).(species).q_rad_el=efluxw_rad_el*n0*te_s0*cs0*ElectronCharge;
  end
end
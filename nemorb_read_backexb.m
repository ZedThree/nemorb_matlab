function [nemo_t]=nemorb_read_backexb(nemo_t, ind)
%######################################
%[nemo_t]=nemorb_read_backexb(nemo_t, info, ind)
%######################################
%---------------
%
%Description
%
% Read background rotation data
%
%---------------
%
%Input arguments
%
% nemo_t = array containing simulation data structure
% ind    = array containing indices of simulations to read equil data from
%---------------
%
%Output
%
% nemo_t = simulation data structure with background rotation data
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

  disp('reading background rotation data')
  vexb_chi  = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/vexb_chi');
  vexb_phi  = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/vexb_phi'); 
  v_par     = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/v_par');    
  omega     = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/omega_prof');
  v_par_pol = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/v_par_pol'); 
  v_par_tor = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/v_par_tor'); 
  Er	    = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/Er'); 
  potential = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/potential');
  time      = hdf5read(nemo_t(k).filename, '/equil/profiles/generic/back_exb/time');   

  nemo_t(k).generic.vexb_chi  = vexb_chi ;
  nemo_t(k).generic.vexb_phi  = vexb_phi ;
  nemo_t(k).generic.v_par     = v_par    ;
  nemo_t(k).generic.omega_prof= omega    ;
  nemo_t(k).generic.v_par_pol = v_par_pol;
  nemo_t(k).generic.v_par_tor = v_par_tor;
  nemo_t(k).generic.Er	      = Er       ;
  nemo_t(k).generic.potential = potential;
  nemo_t(k).generic.time      = time     ;

end % Loop on sims
nemo_t=nemo_t;
cd(pwd_old);

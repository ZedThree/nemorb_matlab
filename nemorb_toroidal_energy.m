function [sim] = nemorb_toroidal_energy(sim,	info, species, graphs, ind)
%################################
%[sim gamma_graph gamma_avg gamma_std] = nemorb_growth_convergance(sim,info, species, graphs, ind)
%################################
%---------------
%
% Description
%
% Plot growth rate of toroidal mode vs. time (intended for single modes)
%
%---------------
%
% Input arguments
%
% sim     = array containing the simulations data, created with nemorb_load.m
% info    = array containing info from nemorb_load
% species = name of species to plot growth rates for
% graphs  = [ 1 : 0 ] graphs on or off
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
% Output
%
% sim         
% gamma_graph = growth from linear fit to log plot of energy
% gamma_avg   = growth rate from average of growth rate evolution
%---------------

global timelab
if exist('ind') == 0
  ind=1;
end

warning('off','curvefit:fit:nonDoubleXData');

pwd_old=pwd;
for i=1:length(ind) %Loop on simulations
  k=ind(i);
  cd(sim(k).path);

  % For backwards compatibility - check where background_efield_size
  % is in the hdf5 file, then make a local variable.
  background_efield_size = nemorb_genericornot(sim(k),info(k),'background_efield_size');
  nfilt1 = nemorb_genericornot(sim(k),info(k),'nfilt1');
  lx = nemorb_genericornot(sim(k),info(k),'lx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thermal velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ion_mass = 2*(1.6726E-27);
  charge = 1.6022E-19;
  % k_B = 1.3807E-23;
  a = sim(k).generic.a_mid/sim(k).generic.d_norm;
  if sim(k).(species).nsel_profile == 5
    speak_ind=find(sim(k).(species).s_prof>sim(k).(species).speak,1);
    kappan = sim(k).(species).gradn_pic(speak_ind);
  else
    kappan = sim(k).(species).kappan0;
  end
  L_n = kappan/a;
  omega_ci = charge*sim(k).generic.btor0/ion_mass;
  mass = sim(k).(species).mass;
  gamma_norm = (lx/2)*sqrt(mass);%/abs(kappan);%*sqrt(1)));%/abs(kappan);
  %gamma_norm = 1;
				 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check dataset isn't empty - matlab will crash if attempts to load
%an empty dataset
  efspec=strcat('/data/var2d/',species,'/efspec1mn');
  %efpsec=['/data/vard2d/generic/phispecmn']
  efspec_time=strcat(efspec,'/time');
  dims=hdf5_checkdims(info,efspec_time,ind);
  if dims==0
    err=sprintf('%s','ERROR: efspec1mn is empty, quitting');
    disp(err)
    return
  end
  %Read poloidal and toroidal modes
  efspec_coord1=strcat(efspec,'/coord1');
  efspec_coord2=strcat(efspec,'/coord2');
  n=hdf5read(sim(k).filename,efspec_coord1);
  m=hdf5read(sim(k).filename,efspec_coord2);
  nmodes=length(n);
  mmodes=length(m);

  time_ef=hdf5read(sim(k).filename,efspec_time);
  nsteps_ef=length(time_ef);
  efspecn=zeros(nmodes,nsteps_ef);
  dt_ef=time_ef(2)-time_ef(1);


  for j=1:length(n) %Loop on toroidal modes
    nloc=n(j);
    clear tmp;
    %Get efspec(nloc,:,:)= (m,t) components of toroidal mode nloc
    %warning in hdf5read_slice: dimensions in reverse order!!!
    efspec_data=strcat(efspec,'/data');
    tmp=hdf5read_slice_new(sim(k).filename,efspec_data,[0 0 j-1],[nsteps_ef mmodes 1]);
    %remove the 1st useless dimension
    tmp=squeeze(tmp);
    %tmp is a 2d array (m,t)
    %Sum over poloidal modes to get the total toroidal energy
    efspecn(j,:)=abs(sum(tmp,1));
    leg{j}=strcat('n=',num2str(nloc));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end %Loop on toroidal modes

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (graphs == 1)
    % Plot energy of toroidal modes
    plotname = [sim(k).name ' - Energy'];
    figure('Name',plotname);
    h=surf(time_ef,n,log10(efspecn));
    shading interp;
%    set(get(h,'Parent'),'ZScale','log');
    xlabel(timelab{1})
    ylabel('E^{(n)} [m_ic_s^2]')
    title(sprintf('Energy  %s',sim(k).name))
  end

end %loop on simulations

sim=sim;
cd(pwd_old);
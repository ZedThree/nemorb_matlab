function [sim gamma_graph gamma_avg gamma_std] = nemorb_growth_convergance(sim, ...
				info, species, graphs, ind)
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

  %Check that linear times for growth rate calculation are inside time interval
  if sim(k).tmax > sim(k).time(end)
    disp('tmax modified because it was too big')
    sim(k).tmax=sim(k).time(end);
  end

  if sim(k).tmin > sim(k).time(end)
    disp('tmin modified because it was too big')
    sim(k).tmin=sim(k).time(1);
  end
  imin=1+floor(sim(k).tmin/dt_ef);
  imax =1+floor(sim(k).tmax/dt_ef);

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
    Estart(j)=(efspecn(j,imin));
    Estop(j) =(efspecn(j,imax));

    %%%%%%%% test %%%%%%%%
    %warning('off','optim:lsqlin:LinConstraints')
    [a c]=fit(time_ef(imin:imax),log(efspecn(j,imin:imax))','poly1');
    sim(k).fit_coefs=a;
    sim(k).fit_goodness=c;
    xx=imin;xy=imax;
    test=coeffvalues(a);
    sim(k).gamma2(j)=test(1)*gamma_norm/2;
    ci=confint(a);
    sim(k).gamma_confint=ci(:,1)*gamma_norm/2;
    test_values=test(1)*time_ef(imin:imax)+test(2);
    test_time=time_ef(imin:imax);
    residual=efspecn(j,imin:imax)-exp(test_values');
    %test_values=test(2)*exp(test(1)*time_ef(imin:imax));
    %%%%%%%% test %%%%%%%%

    %Get linear growth rate
    [coef, goodness]=polyfit(double(time_ef(imin:imax)),transpose(log(efspecn(j,imin:imax))),1);
    sim(k).gamma(j)=gamma_norm*coef(1)/2;
    %Get error on gamma
    ste = sqrt(diag(inv(goodness.R)*inv(goodness.R')).*goodness.normr.^2./goodness.df); 
    %3D sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df)
    %[y delta]=polyval(coef,double(time_ef(imin:imax)),fit);
    sim(k).gamma_ste(j)=ste(1)*gamma_norm/2;

%%%%%%%%%% Another way - fitting an exp %%%%%%%%%%%%
% $$$ X0 = [1 1 0.01];
% $$$ options = optimset('Largescale','off');
% $$$ x=lsqnonlin(@exp_fit,X0,[],[],options,double(time_ef(imin:imax)),efspecn(j,imin:imax));
% $$$ Y_new(k,:)=x(1) + x(2).*exp(x(3).*time_ef(imin:imax));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end %Loop on toroidal modes

  %Fastest growing mode
  [gmax n0]=max(sim(k).gamma);
  n0max=n(n0);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Growth rate evolution  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Evolution of growth rate - take the growth rate of small steps -
  %for linear modes, the growth rate should converge
  jj = 1;
  gamma_cov = 0.0;			% convergant growth rate
  time_growth = 0.0;
  time_step=100;
  %for ii = 0:dt_ef:(time_ef(end)-time_step)
  for ii = 0:dt_ef:(sim(k).tmax-time_step)
    temp_tmin = ii;
    temp_tmax = ii + time_step;
    imin = 1 + floor(temp_tmin/dt_ef);
    imax = 1 + floor(temp_tmax/dt_ef);
    [coef, goodness]=polyfit(double(time_ef(imin:imax)), ...
			     transpose(log(efspecn(1,imin:imax))),1);
    gamma_cov(jj)=gamma_norm*coef(1)/2;
    time_i = 1 + floor(imin+(imax-imin)/2);
    time_growth(jj)=time_ef(time_i);
    jj=jj+1;
  end
  % average gamma_cov from sim.tmin to sim.tmax
  time_avg_min = 1 + floor(sim(k).tmin/dt_ef);
  jj_avg_min = 1 + floor(time_growth(time_avg_min)/dt_ef);
  time_avg_max = floor((sim(k).tmax-time_step)/dt_ef) - 1;

  if time_avg_max < 1, time_avg_max = size(time_growth,2); end
  jj_avg_max = floor(time_growth(time_avg_max)/dt_ef) - 1;
  if jj_avg_max > size(gamma_cov,2), jj_avg_max = size(gamma_cov,2); end

  sim(k).gamma_avg = mean(gamma_cov(jj_avg_min:jj_avg_max));
  sim(k).gamma_std = std(gamma_cov(jj_avg_min:jj_avg_max));
  gamma_graph(k,:) = sim(k).gamma;
  gamma_avg(k,:) = sim(k).gamma_avg;
  gamma_std(k,:) = sim(k).gamma_std;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (graphs == 1)
    % Plot energy of toroidal modes
    plotname = [sim(k).name ' - Growth rate'];
    figure('Name',plotname);
    if nfilt1 == 0
      semilogy(time_ef,(efspecn(1,:)),'k--');
      hold on;
      choose_n = [ 20 40 60 80];
      semilogy(time_ef,(efspecn(choose_n,:)))
%       semilogy(time_ef,(efspecn(2:end,:)))
    else
      semilogy(time_ef,efspecn(1:end,:))
      hold on;
    end
    xlabel(timelab{1})
    ylabel('E^{(n)} [m_ic_s^2]')
    title(sprintf('Energy in mode n=%d, gamma=%5.3f, for %s',n(n0),sim(k).gamma(n0),sim(k).name))
    %legend(leg,'Location','NorthEast')
    
%     semilogy(time_ef(xx:xy),exp(test_values),'g--')
    
%     t=[sim(k).tmin sim(k).tmax];
%     E=[Estart; Estop];
%     plot(t,E,'r:x');
    
    lplot2=0;
    if lplot2==1
      t_growth=[time_growth(time_avg_min), time_growth(time_avg_max)];
      gamma_temp=[gamma_avg(k), gamma_avg(k)];
      
      % Plot evolution of growth rate
      figure;
      hold on;
      plot(time_growth, gamma_cov)
      plot(t_growth, gamma_temp ,'r--')
      axis ([time_ef(1), time_ef(end),  min(gamma_cov)-0.005, max(gamma_cov)+0.005]);
      xlabel(timelab{1})
      %ylabel('Growth rate [v_{th}/L_n]')
      ylabel('Growth rate [v_{th}/a]')
      title(sprintf('Growth rate evolution of %s with gamma=%0.3g',sim(k).name,sim(k).gamma(n0) ))
    end
  end

  sprintf('%s has gamma = %d, and averaged gamma = %d' ...
	  , sim(k).name,sim(k).gamma, sim(k).gamma_avg);

  %Plot growth rates of toroidal modes
% $$$ figure;
% $$$ plot(n,sim(k).gamma,'bx-')
% $$$ xlabel('n')
% $$$ ylabel('\gamma [\Omega_i]')
end %loop on simulations

sim=sim;
cd(pwd_old);
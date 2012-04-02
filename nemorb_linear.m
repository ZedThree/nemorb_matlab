function [gamma sigma_g gamma_inst]=nemorb_linear(nemo,s0d,lplot,ind)
%######################################
%[gamma sigma_g]=nemorb_linear(nemo,s0d,lplot,ind)
%######################################
%---------------
%
%Description
%
%For linear runs, compute growth rate of toroidal mode
%
%---------------
%
%Dummy arguments
%
%nemo   = array containing the simulations data, created with loader.m
%s0d   = species label
%lplot : true = plot results
%ind   = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%gamma = growth rate in \Omega_i^{-1} units
%error on growthrate
%---------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GROWTH RATE (LINEAR SIMULATIONS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  global timelab

  if exist('nemo')==0
    help nemorb_linear;
  else

  if exist('ind')==0
    ind=1;
  end

  pwd_old=pwd;
  colorp={'blue','red','black','cyan','yellow','green'};
  for ii=1:length(ind) %loop on simulations
    k=ind(ii);
    %go to simulation directory
    cd(nemo(k).path);
    %
    ninst=hdf5read(nemo(k).filename,['/parameters/generic/nfreq_spec']);
    ninst=5
    %read n and m modes
    n=hdf5read(nemo(k).filename,['/data/var2d/' s0d '/efspec1mn/coord1']);
    m=hdf5read(nemo(k).filename,['/data/var2d/' s0d '/efspec1mn/coord2']);
    nmodes=length(n);
    mmodes=length(m);
    %read associated time
    tmp =hdf5read(nemo(k).filename,['/data/var2d/' s0d '/efspec1mn/time']);
    nsteps_ef(ii)=length(tmp);
    time_ef(1:nsteps_ef(ii),ii)=tmp;
    dt_ef=time_ef(2,ii)-time_ef(1,ii);
    %read time for phispec
    tmp=hdf5read(nemo(k).filename,['/data/var2d/generic/phispecmn/time']);
    nsteps_phi(ii)=length(tmp);
    time_phi(1:nsteps_phi(ii),ii)=tmp;
    dt_phi=time_phi(2,ii)-time_phi(1,ii);
    mass=nemo(k).(s0d).mass;
    %check that the boundary times to compute the growth rates are correct
    %
    if nemo(k).tmax > nemo(k).time(end)
      nemo(k).tmax=nemo(k).time(end);
      disp(sprintf('tmax modified because it was too big, new tmax = %15.5f', nemo(k).tmax))
    end
    if nemo(k).tmin > nemo(k).time(end)
      nemo(k).tmin=nemo(k).time(1);
      disp(sprintf('tmin modified because it was too big, new tmin = %15.5f', nemo(k).tmin))
    end
    imin_ef=1+floor(nemo(k).tmin/dt_ef);
    imax_ef=1+floor(nemo(k).tmax/dt_ef);
    imin_phi=1+floor(nemo(k).tmin/dt_phi);
    imax_phi=1+floor(nemo(k).tmax/dt_phi);
    clear tmp;
    %Get efspec(nloc,:,:) 
    %warning in hdf5read_slice: dimensions in reverse order!!!
    %Read E(n0,m,t)
    tmp=hdf5read_slice(nemo(k).filename,['/data/var2d/' s0d '/efspec1mn/data'],[0 0 0],[nsteps_ef(ii) mmodes 1]);
    tmp2=hdf5read_slice(nemo(k).filename,['/data/var2d/generic/phispecmn/data'],[0 0 0],[nsteps_phi(ii) mmodes 1]);
    %tmp=tmp%+tmp3;
    %remove the 1st useless dimension (n)
    tmp=squeeze(tmp);
    tmp2=squeeze(tmp2);
    %sum over poloidal modes to get the total toroidal mode energy
    tmp=squeeze(sum(tmp,1));
    tmp2=squeeze(sum(tmp2,1));
    %
    efspecn(1:nsteps_ef(ii),ii)=abs(tmp);
    efspecnlog(1:nsteps_ef(ii),ii)=log(efspecn(1:nsteps_ef(ii),ii));
    %
    phispecn(1:nsteps_phi(ii),ii)=tmp2;
    phispecnlog(1:nsteps_phi(ii),ii)=log(phispecn(1:nsteps_phi(ii),ii));
    leg{ii}=strcat('n=',num2str(j-1),', ',nemo(k).name);
    %extract growth rate with a linear fit
    [coef, fit]=polyfit(double(time_ef(imin_ef:imax_ef,ii)),efspecnlog(imin_ef:imax_ef,ii),1);
    [coef2, fit2]=polyfit(double(time_phi(imin_phi:imax_phi,ii)),phispecnlog(imin_phi:imax_phi,ii),1);
    %divide by 2: field energy ~ exp(2*gamma*t)
    gamma(ii)=coef(1)/2;
    gamma2(ii)=coef2(1)/2;
    t(1,ii)=nemo(k).tmin; 
    t(2,ii)=nemo(k).tmax;

    %instantaneous growthrate

    out=0;
    tt=1;
    itt1=1;
    itt2=1;
    while out==0
      itt1=itt2;
      itt2=itt1+ninst;
      if itt2 < nsteps_ef(ii)
	gamma_inst(tt,ii)=(efspecnlog(itt2,ii)-efspecnlog(itt1,ii))/(time_ef(itt2,ii)-time_ef(itt1,ii));
	time_inst(tt,ii)=(time_ef(itt2,ii)+time_ef(itt1,ii))/2.0;
	tt=tt+1;
      else
	out=1;
      end
    end

    ntime_inst(ii)=tt-1;
  end %loop on simulation

  gamma_inst=gamma_inst*0.5;

  %get error on growthrate
  for ii=1:length(ind)
    k=ind(ii);
    its=find(time_inst(:,ii) >=nemo(k).tmin,1,'first');
    ite=find(time_inst(:,ii) >=nemo(k).tmax,1,'first');
    if isempty(ite) ==1
      ite=size(time_inst,1);
    end
    if time_inst(end,ii) == 0
      ite=find(time_inst(:,ii) ~= 0.0,1,'last');
    end


    sigma_g(ii)=std(gamma_inst(its:ite,ii));
    %global growth rate = average of instantaneous growth rate
    mean_gamma(ii)=mean(gamma_inst(its:ite,ii));
    disp(sprintf('global growthrate E_f= %0.5e',gamma(ii)))
    disp(sprintf('global growthrate (phi^2) = %0.5e',gamma2(ii)))
    disp(sprintf('mean of instantaneus growthrate = %0.5e',mean_gamma(ii)))
    disp(sprintf('error = %0.5e',sigma_g(ii)))
  end
  if lplot
    %plot result
    figure;
    for ii=1:length(ind)
      k=ind(ii);
      ic=mod(ii,length(colorp));
      leg{ii}=strcat('n=',num2str(nemo(k).generic.nfilt1),', sim=',nemo(k).name);
      semilogy(time_phi(1:nsteps_phi(ii),ii),phispecn(1:nsteps_phi(ii),ii),'-','Color',colorp{ic});
      hold on
      semilogy(time_phi(imin_phi:imax_phi,ii), exp(coef2(2)+coef2(1)*time_phi(imin_phi:imax_phi,ii)),'-','Color',colorp{ic+1});
      semilogy(time_phi(imin_phi:imax_phi,ii), exp(coef2(2)+2*mean_gamma(ii)*time_phi(imin_phi:imax_phi,ii)),'-','Color',colorp{ic+2});
      xlabel(timelab{1})
      ylabel(strcat('phi^{(n)}'))
    end
    semilogy([nemo(k).tmin nemo(k).tmin], [min(log10(phispecn(:,ii))) max(log10(phispecn(:,ii)))],'-','Color','black');
    semilogy([nemo(k).tmax nemo(k).tmax], [min(log10(phispecn(:,ii))) max(log10(phispecn(:,ii)))],'-','Color','black');
    legend([leg 'fit' 'instantaneus'], 'Location','NorthWest')
    figure;
    for ii=1:length(ind)
      k=ind(ii);
      ic=mod(ii,length(colorp));
      leg{ii}=strcat('n=',num2str(nemo(k).generic.nfilt1),', sim=',nemo(k).name);
      semilogy(time_ef(1:nsteps_ef(ii),ii),efspecn(1:nsteps_ef(ii),ii),'-','Color',colorp{ic});
      hold on
      semilogy(time_ef(imin_ef:imax_ef,ii), exp(coef(2)+coef(1)*time_ef(imin_ef:imax_ef,ii)),'-','Color',colorp{ic+1});
      semilogy(time_ef(imin_ef:imax_ef,ii), exp(coef(2)+2*mean_gamma(ii)*time_ef(imin_ef:imax_ef,ii)),'-','Color',colorp{ic+2});
      xlabel(timelab{1})
      ylabel(strcat('E^{(n)} [m_ic_s^2]'))
    end
    semilogy([nemo(k).tmin nemo(k).tmin], [min(log10(efspecn(:,ii))) max(log10(efspecn(:,ii)))],'-','Color','black');
    semilogy([nemo(k).tmax nemo(k).tmax], [min(log10(efspecn(:,ii))) max(log10(efspecn(:,ii)))],'-','Color','black');
    legend([leg 'fit' 'instantaneus'], 'Location','NorthWest')
    figure;
     for ii=1:length(ind)
       k=ind(ii);
       ic=mod(ii,length(colorp))+1;
       leg{ii}=strcat('n=',num2str(nemo(k).generic.nfilt1),', sim=',nemo(k).name);
       plot(time_inst(1:ntime_inst(ii),ii),gamma_inst(1:ntime_inst(ii),ii)*mass,'-','Color',colorp{ic});
       xlabel(timelab{1})
       ylabel('\gamma [\Omega_i]')
       hold on
     end
     for ii=1:length(ind)
       ic=mod(ii,length(colorp))+1;
       plot([nemo(k).tmin nemo(k).tmax],[gamma(ii) gamma(ii)],'--','Color',colorp{ic+1})
       errorbar([nemo(k).tmin nemo(k).tmax],[mean_gamma(ii) mean_gamma(ii)],[sigma_g(ii) sigma_g(ii)], '--','Color',colorp{ic+2})
     end
     legend([leg 'fit' 'instantaneus'], 'Location','NorthWest')
  end

  cd(pwd_old);
end

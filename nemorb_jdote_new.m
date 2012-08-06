function [sim1_new,gamma_je]=nemorb_jdote_new(sim1,ind, tmin, tmax)
%######################################
%sim1_new=jdote_hdf5(sim1,s0d,ind)
%######################################
%---------------
%
%Description
%
%For linear runs:
%compute and plot components of energy transfer
%obtain growth rate via energy transfer
%plot power balance
%
%---------------
%
%Dummy arguments
%
%sim1  = array containing the simulations data, created with nemorb.m
%ind   = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
% gamma_je: growth rate estimate from J.dot.E diagnostics
%---------------

  global timelab

  pwd_old=pwd;
  if exist('ind')==0
    ind=1;
  end

  nsimul=length(ind); % Number of simulations

  for i=1:nsimul %--- Loop on simulations ---
    k=ind(i);
    n=sim1(k).nsteps % Number of timesteps of the simulation
    nrspec=length(sim1(k).speciesname)-1 % Nunmber of species, (NB: -1 for 'generic')
    jdote_tot_all=zeros(n,nsimul);
    efield_all=zeros(n,nsimul);

    for jspec=1:nrspec %--- Loop on species ---

      specname=sim1(k).speciesname(jspec);
      data={'jdote_tot','jdote_par','jdote_grb','jdote_curv','jdote_grp','jdote_ebg', 'jdote_exb', 'efield','etransw'};
      sim1=read_0d_hdf5(sim1,specname,data,ind);
      if isfield(sim1(ind).speciesname(jspec),'jdote_exb')==1
        data={'jdote_exb'};
        sim1=read_0d_hdf5(sim1,specname,data,ind);
      else
        sim1.specname.jdote_exb(1:n)=0.;
      end

% fill in jdote(1:n,i) components
      jdote_tot(1:n,jspec,i) = sim1(k).jdote_tot(1:n);
      jdote_par(1:n,jspec,i) = sim1(k).jdote_par(1:n);
      jdote_grb(1:n,jspec,i) = sim1(k).jdote_grb(1:n);
      jdote_curv(1:n,jspec,i) = sim1(k).jdote_curv(1:n);
      jdote_grp(1:n,jspec,i) = sim1(k).jdote_grp(1:n);
      jdote_ebg(1:n,jspec,i) = sim1(k).jdote_ebg(1:n);
      jdote_exb(1:n,jspec,i) = sim1(k).jdote_exb(1:n);
% fill in field energy
      efield(1:n,jspec,i) = sim1(k).efield(1:n);
% sum over species
      jdote_tot_all(1:n,i)=jdote_tot_all(1:n,i)+jdote_tot(1:n,jspec,i);
      efield_all(1:n,i)=efield_all(1:n,i)+efield(1:n,jspec,i);

    end %--- Loop on species ---

% Compute instantaneous growth rates

% (1/2E_field) (dE_field / dt) using finite differences at mid intervals
    dEfdt(1:n-1,i) = (efield_all(2:n)-efield_all(1:n-1))/sim1(k).generic.dt;
    efield_all_mid(1:n-1,i) = 0.5*(efield_all(2:n)+efield_all(1:n-1));
    time_mid(1:n-1,i) = 0.5*(sim1(k).time(2:n)+sim1(k).time(1:n-1));
    for jspec=1:nrspec
      jdote_tot_mid(1:n-1,jspec,i) = 0.5*(jdote_tot(2:n,jspec,i)+jdote_tot(1:n-1,jspec,i));
      jdote_par_mid(1:n-1,jspec,i) = 0.5*(jdote_par(2:n,jspec,i)+jdote_par(1:n-1,jspec,i));
      jdote_grb_mid(1:n-1,jspec,i) = 0.5*(jdote_grb(2:n,jspec,i)+jdote_grb(1:n-1,jspec,i));
      jdote_curv_mid(1:n-1,jspec,i) = 0.5*(jdote_curv(2:n,jspec,i)+jdote_curv(1:n-1,jspec,i));
      jdote_grp_mid(1:n-1,jspec,i) = 0.5*(jdote_grp(2:n,jspec,i)+jdote_grp(1:n-1,jspec,i));
      jdote_ebg_mid(1:n-1,jspec,i) = 0.5*(jdote_ebg(2:n,jspec,i)+jdote_ebg(1:n-1,jspec,i));
      jdote_exb_mid(1:n-1,jspec,i) = 0.5*(jdote_exb(2:n,jspec,i)+jdote_exb(1:n-1,jspec,i));
    end
    gamma_efield = 0.5 * dEfdt ./ efield_all_mid;
    gamma_tot_all=zeros(n-1,nsimul);
    for jspec=1:nrspec
      gamma_tot(1:n-1,jspec,i) = 0.5 * (-jdote_tot_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_par(1:n-1,jspec,i) = 0.5 * (-jdote_par_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_grb(1:n-1,jspec,i) = 0.5 * (-jdote_grb_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_curv(1:n-1,jspec,i) = 0.5 * (-jdote_curv_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_grp(1:n-1,jspec,i) = 0.5 * (-jdote_grp_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_ebg(1:n-1,jspec,i) = 0.5 * (-jdote_ebg_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_exb(1:n-1,jspec,i) = 0.5 * (-jdote_exb_mid(1:n-1,jspec,i)) ./ efield_all_mid(1:n-1,i);
      gamma_tot_all(1:n-1,i) = gamma_tot_all(1:n-1,i) + gamma_tot(1:n-1,jspec,i);
    end

% Determine time interval for fits and averages
    if exist('tmin')==0 % Use hdf5 data for tmin if not given as argument
      tmin=sim1(k).tmin;
    end
    if exist('tmax')==0 % Use hdf5 data for tmax if not given as argument
      tmax=sim1(k).tmax;
    end
% Check bounds
    tmin=min(tmin,sim1(k).time(end))
    tmax=min(tmax,sim1(k).time(end))
% Indices of time interval for fits and averages
    imin=1+floor(tmin/sim1(k).generic.dt);
    imax=1+floor(tmax/sim1(k).generic.dt);

% Fit -J.E(t) ~ j0*exp(gamma*t)
    [coefj, fitj]=polyfit(double(sim1(k).time(imin:imax)),double(log(abs(jdote_tot_all(imin:imax,i)))),1);
    gamma_jdote(i)=coefj(1)/2;

% Fit Efield(t) ~ E0*exp(gamma*t)
    [coefe, fite]=polyfit(double(sim1(k).time(imin:imax)),double(log(abs(efield_all(imin:imax,i)))),1);
    gamma_efield(i)=coefe(1)/2;

% Check power balance growth rate, i.e. -jdote/2E should be equal to gamma (on fits)
    gamma_powbal_fit_tmin=0.5*exp(coefj(2)+coefj(1)*tmin)/exp(coefe(2)+coefe(1)*tmin);
    gamma_powbal_fit_tmax=0.5*exp(coefj(2)+coefj(1)*tmax)/exp(coefe(2)+coefe(1)*tmax);

% Averages of instantaneous growth rates
    gamma_efield_av(i)= mean(gamma_efield(imin:imax-1,i));
    gamma_je_av(i)= mean(gamma_tot_all(imin:imax-1,i)); % -1 because array defined at mid intervals

% Units conversion factors
    specname=sim1(k).speciesname(1); % !!! We use density gradient length of 1st species to obtain Ln !!!
    spname=specname{1};
    Ln_cs=0.5*sim1(k).generic.lx/sim1(k).(spname).kappan0;
    a_cs=0.5*sim1(k).generic.lx;

% Display
    disp( [sprintf('gamma from Ef fit t in [%0.5g',tmin) ',' sprintf('%0.5g]',tmax) sprintf(' [Omega_i^-1]   = %0.8g', gamma_efield(i)), ' [Omega_i]' ] )
    disp( [sprintf('gamma from J.E fit t in [%0.5g',tmin) ',' sprintf('%0.5g]',tmax) sprintf(' [Omega_i^-1]   = %0.8g', gamma_jdote(i)), ' [Omega_i]' ] )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmin) sprintf(' [Omega_i^-1]   = %0.8g', gamma_powbal_fit_tmin), ' [Omega_i]' ]   )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmax) sprintf(' [Omega_i^-1]   = %0.8g', gamma_powbal_fit_tmax), ' [Omega_i]' ]   )
    disp( [sprintf('average gamma from Ef t in [%0.5g',tmin) ',' sprintf('%0.5g]',tmax) sprintf(' [Omega_i^-1]   = %0.8g', gamma_efield_av(i)), ' [Omega_i]' ] )
    disp( [sprintf('average gamma from J.E t in [%0.5g',tmin) ',' sprintf('%0.5g]',tmax) sprintf(' [Omega_i^-1]   = %0.8g', gamma_je_av(i)), ' [Omega_i]' ] )
    disp([' ']) %---------------------------------------
    disp( [sprintf('gamma from Ef fit t in [%0.5g',tmin/Ln_cs) ',' sprintf('%0.5g]',tmax/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_efield(i)*Ln_cs), ' [cs/Ln]' ] )
    disp( [sprintf('gamma from J.E fit t in [%0.5g',tmin/Ln_cs) ',' sprintf('%0.5g]',tmax/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_jdote(i)*Ln_cs), ' [cs/Ln]' ] )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmin/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_powbal_fit_tmin*Ln_cs), ' [cs/Ln]' ]   )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmax/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_powbal_fit_tmax*Ln_cs), ' [cs/Ln]' ]   )
    disp( [sprintf('average gamma from Ef t in [%0.5g',tmin/Ln_cs) ',' sprintf('%0.5g]',tmax/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_efield_av(i)*Ln_cs), ' [cs/Ln]' ] )
    disp( [sprintf('average gamma from J.E t in [%0.5g',tmin/Ln_cs) ',' sprintf('%0.5g]',tmax/Ln_cs) sprintf(' [Ln/cs]   = %0.8g', gamma_je_av(i)*Ln_cs), ' [cs/Ln]' ] )
    disp([' ']) %----------------------------------------
    disp( [sprintf('gamma from Ef fit t in [%0.5g',tmin/a_cs) ',' sprintf('%0.5g]',tmax/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_efield(i)*a_cs), ' [cs/a]' ] )
    disp( [sprintf('gamma from J.E fit t in [%0.5g',tmin/a_cs) ',' sprintf('%0.5g]',tmax/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_jdote(i)*a_cs), ' [cs/a]' ] )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmin/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_powbal_fit_tmin*a_cs), ' [cs/a]' ]   )
    disp( [sprintf('J.E/2Ef fit at t= %0.5g',tmax/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_powbal_fit_tmax*a_cs), ' [cs/a]' ]   )
    disp( [sprintf('average gamma from Ef t in [%0.5g',tmin/a_cs) ',' sprintf('%0.5g]',tmax/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_efield_av(i)*a_cs), ' [cs/a]' ] )
    disp( [sprintf('average gamma from J.E t in [%0.5g',tmin/a_cs) ',' sprintf('%0.5g]',tmax/a_cs) sprintf(' [a/cs]   = %0.8g', gamma_je_av(i)*a_cs), ' [cs/a]' ] )

% Plot field energy vs time
    figure
    hold on
    h=plot(sim1(k).time, efield_all(1:n,i),'b-');
    xlabel(timelab{1})
    ylabel('E_f')
    set(gca,'yscale','log')
    % superpose fit
    hfit=plot(sim1(k).time, exp(coefe(2)+coefe(1)*sim1(k).time),'b--');
    title(['Field energy ', sim1(k).name])
    box on
    % superpose fit
    hfit=plot(sim1(k).time, exp(coefe(2)+coefe(1)*sim1(k).time),'b--');
    legend('E_f','fit','Location','NorthWest')
 
% Plot -J.E vs time
    figure
    hold on
    h=plot(sim1(k).time,+jdote_tot_all(1:n,i),'r--');
    h=plot(sim1(k).time,-jdote_tot_all(1:n,i),'k-');
    xlabel(timelab{1})
    ylabel('|J.E|')
    set(gca,'yscale','log')
    specname=sim1(k).speciesname(nrspec+1); spname=specname{1};
    strtitle=['Field-particle power ', sim1(k).name,' ', spname];
    title(strtitle);
    box on
     % superpose fit
    hfit=plot(sim1(k).time, exp(coefj(2)+coefj(1)*sim1(k).time),'k--');
    legend('+J.E','-J.E','fit','Location','NorthWest');
   if (nrspec>1)
      for jspec=1:nrspec
	figure
	hold on
	h=plot(sim1(k).time,+jdote_tot(1:n,jspec,i),'r--');
	h=plot(sim1(k).time,-jdote_tot(1:n,jspec,i),'k-');
	xlabel(timelab{1})
	ylabel('|J.E|')
	set(gca,'yscale','log')
	specname=sim1(k).speciesname(jspec); spname=specname{1};
	strtitle=['Field-particle power ', sim1(k).name,' ', spname];
	title(strtitle)
	box on
	legend('+J.E','-J.E','Location','NorthWest')
       end
    end

% Plot separate contributions to instantaneous growth rate
  for jspec=1:nrspec
    figure
    hold on
    h=plot(time_mid(1:n-1,i),gamma_tot(1:n-1,jspec,i),'k-');
    xlabel(timelab{1})
    ylabel('\gamma [\Omega_i]')
    h=plot(time_mid(1:n-1,i),gamma_par(1:n-1,jspec,i),'b-');
    h=plot(time_mid(1:n-1,i),gamma_grb(1:n-1,jspec,i),'r-');
    h=plot(time_mid(1:n-1,i),gamma_curv(1:n-1,jspec,i),'c-');
    h=plot(time_mid(1:n-1,i),gamma_grp(1:n-1,jspec,i),'g-');
    h=plot(time_mid(1:n-1,i),gamma_exb(1:n-1,jspec,i),'m-');
    h=plot(time_mid(1:n-1,i),gamma_ebg(1:n-1,jspec,i),'m--');
    specname=sim1(k).speciesname(jspec); spname=specname{1};
    strtitle=['Contributions to \gamma ', sim1(k).name,' ', spname];
    title(strtitle)
    box on
    legend('tot','//','grad-B','curv','grad-P','ExB','E_bxB','Location','NorthWest')
  end

% Plot instantaneous growth rates from Efield and J.E
  figure
  hold on
  h=plot(time_mid(1:n-1,i),gamma_tot_all(1:n-1,i),'k-');
  h=plot(time_mid(1:n-1,i),gamma_efield(1:n-1,i),'b-');
  xlabel(timelab{1})
  ylabel('\gamma [\Omega_i]')
  strtitle=['\gamma from Efield and from J.E ',sim1(k).name];
  title(strtitle)
  box on


  end %--- Loop on simulations



gamma_je=1.2345;
sim1_new=sim1;

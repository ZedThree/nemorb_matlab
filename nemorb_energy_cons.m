function nemo_new=nemorb_energy_cons(nemo,s0d,ind)
%################################
%nemo_new=energy_cons_hdf5(nemo,s0d,ind)
%################################
%---------------
%
%Description
%
%Plot the relative and absolute conservation of energy
%
%---------------
%
%Dummy arguments
%
%nemo = array containing the simulations data, created with nemorb.m
%s0d = species name
%ind = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%sim_new = simulation data with added conservation energy diagnostic
%---------------
%
  global timelab
  pwd_old=pwd;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENERGY CONSERVATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('nemo')==0
    help nemorb_energy_cons;
  else
    if exist('ind')==0
      ind=1;
    end
    colorp={'blue','red','green','black','cyan','yellow'};
    data={'ekin','ekinf0','efield','etransw','etransf0'};
%    nemo=read_0d_hdf5(nemo,s0d,data,ind);
    nlelec=false;
    if isfield(nemo.electrons, 'efluxf0_rad')==1
      for i=1:length(ind)
	k=ind(i);
	nlelec = false;
	if (nemo(k).electrons.nsel_push_trap_only > 0)
	  nlelec = true;
	end
      end
    end
    if nlelec
      data={'ekin','ekinf0','efield','etransw','etransf0'};
%      nemo=read_0d_hdf5(nemo,electrons,data,ind);
    end    
    for i=1:length(ind) %loop on simulations
      k=ind(i);
      n=nemo(k).nsteps;
      %Total energy
      %E=E_kin0+E_kin+E_field
      if nlelec
	ekin_tot(1:n,i)  = nemo(k).(s0d).ekin   + nemo(k).electrons.ekin;
	ekinf0_tot(1:n,i)= nemo(k).(s0d).ekinf0 + nemo(k).electrons.ekinf0;
	efield_tot(1:n,i)= nemo(k).(s0d).efield + nemo(k).electrons.efield;
	e(1:n,i)         = ekin_tot(1:n,i)     + efield_tot(1:n,i);
	etransw_tot(1:n,i)  = nemo(k).(s0d).etransw +  nemo(k).electrons.etransw;
	etransf0_tot(1:n,i) = nemo(k).(s0d).etransf0 + nemo(k).electrons.etransf0;
      else
	ekin_tot(1:n,i)  = nemo(k).(s0d).ekin;
	ekinf0_tot(1:n,i)= nemo(k).(s0d).ekinf0;
	efield_tot(1:n,i)= nemo(k).(s0d).efield;
	e(1:n,i)=ekin_tot(1:n,i)+ekinf0_tot(1:n,i)+efield_tot;
	%Total transferred energy
	etransw_tot(1:n,i)=nemo(k).(s0d).etransw;
	etransf0_tot(1:n,i)=nemo(k).(s0d).etransf0;
      end
      %
      %E_f(t)-E_f(t_0)
      delta_efield(1:n,i)=efield_tot-efield_tot(1);
      %E_kin(t)-E_kin(t_0)
      delta_ekin(1:n,i)=ekin_tot(1:n,i)-ekin_tot(1,i);
      %E_kin0(t)-E_kin0(t_0)
      delta_ekinf0(1:n,i)=ekinf0_tot(1:n,i)-ekinf0_tot(1,i);
      %E_kin(t)+E_kin0(t)-E_kin(t_0)-E_Kin0(t_0)
      delta_ekin_tot(1:n,i)=delta_ekin(1:n,i)+delta_ekinf0(1:n,i);
      %
      %Integrate power transfer to get ekin(t)-ekin(t_0)
      dek=zeros(length(nemo(k).(s0d).etransw),length(ind));
      dek0=zeros(length(nemo(k).(s0d).etransf0),length(ind));
      %
      for j=2:nemo(k).nsteps,
	dt_int    = nemo(k).time(j)-nemo(k).time(j-1);
	dek(j,i)  = dek(j-1,i)  + 0.5*dt_int*(etransw_tot(j,i)  + etransw_tot(j-1,i));
	dek0(j,i) = dek0(j-1,i) + 0.5*dt_int*(etransf0_tot(j,i) + etransf0_tot(j-1,i));
      end

      %E_kin(t)+E_kin0(t)-E_kin(t_0)-E_Kin0(t_0)
      delta_ekin_int_tot(1:n,i)=dek(1:n,i)+dek0(1:n,i);

      %dE=E(t)-E(t_0)
      delta_e(1:n,i)=delta_efield(1:n,i)+delta_ekin_tot(1:n,i);
      delta_eint(1:n,i)=delta_efield(1:n,i)+delta_ekin_int_tot(1:n,i);

      %dE/E_f(t)
      relcons(1:n,i)     = delta_e(1:n,i)    ./nemo(k).(s0d).efield;
      relcons_int(1:n,i) = delta_eint(1:n,i) ./nemo(k).(s0d).efield;

      %dE/E(t_0)
      abscons(1:n,i)     = delta_e(1:n,i)     /e(1);
      abscons_int(1:n,i) = delta_eint(1:n,i)  /e(1);
    end

    clear legt leg
    legt{1}='E=E_f+E_k';
    legt{2}='E=E_f+\int E\prime_k dt';

    %Plot relative energy conservation, direct and integrated
    figure;
    hold on
    for i=1:length(ind)
      k=ind(i);
      n=nemo(k).nsteps;
      ic=mod(i,length(colorp))+1;
      plot(nemo(k).time,relcons(1:n,i),'-','Color',colorp{ic})
      plot(nemo(k).time,relcons_int(1:n,i),'--','Color',colorp{ic+1});
      for j=1:2
	m=(i-1)*2+j;
	leg{m}=strcat(legt{j},',',nemo(k).name);
      end
    end
    xlabel(timelab{1})
    ylabel('(E(t)-E(t_0))/E_f(t)')
    legend(leg,'Location','SouthWest')
    ymin=max(max(abs(relcons(end,:)),abs(relcons_int(end,:))));
    axis_old=axis;
    axis([axis_old(1) axis_old(2) -abs(ymin) abs(ymin)])
    title('Relative energy conservation')

    %Plot absolute energy conservation, direct and integrated
    figure;
    hold on;
    for i=1:length(ind)
      k=ind(i);
      n=nemo(k).nsteps;
      ic=mod(i,length(colorp))+1;
      plot(nemo(k).time,abscons(1:n,i),'-','Color',colorp{ic})
      plot(nemo(k).time,abscons_int(1:n,i),'--','Color',colorp{ic+1})
    end
    xlabel(timelab{1})
    ylabel('(E(t)-E(t_0))/E(t_0)')
    legend(leg,'Location','SouthWest')
    title('Absolute energy conservation')

    clear legt leg
    legt{1}='\Delta E_{tot}';
    legt{2}='\Delta E_{tot} with \int E\prime_k';
    legt{3}='\Delta E_k';
    legt{4}='\Delta E_k with \int E\prime_k';
    legt{5}='\Delta E_f';

    %Plot components of energy, direct and integrated
    for i=1:length(ind)
      k=ind(i);
      n=nemo(k).nsteps;
      figure;
      hold on
      plot(nemo(k).time,delta_e(1:n,i),'b-')
      plot(nemo(k).time,delta_eint(1:n,i),'b--')
      plot(nemo(k).time,delta_ekin_tot(1:n,i),'r-')
      plot(nemo(k).time,delta_ekin_int_tot(1:n,i),'r--')
      plot(nemo(k).time,delta_efield(1:n,i),'k')
      legend(legt,'Location','NorthWest')
      xlabel(timelab{1})
      ylabel('\Delta E [m_ic_s^2]')
      title(strcat('Energy deviation , ',nemo(k).name))
    end

    nemo_new=nemo;

    cd(pwd_old)
  end

    
    

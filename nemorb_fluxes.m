function nemo_t=nemorb_fluxes(nemo_t, lplot, tol,ind)
%######################################
%nemo_t=nemorb_fluxes(nemo_t, lplot, tol, ind)
%######################################
%---------------
%
%Description
%
%Compute different transport quantities for ions
%
%---------------
%
%Input arguments
%
%nemo_t   = array containing the nemo_tulations data, created with loader.m
%lplot  : true = plot reconstructed density and temperature
%tol   = tolerance for splines smoothing
%ind   = array containing indexes of nemo_tulations for the output of this function
%---------------
%
%Output
%
%nemo_t = new structure with the 1D quantities loaded in this script
%---------------

global timelab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HEADERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('ind')==0
  ind=1;
end
% Turn off warnings about NaNs in smoothing
warning('off','SPLINES:CHCKXYWP:NaNs');
colorp={'blue','red','green','black','cyan','yellow'};
for ii=1:length(ind)
  kk=ind(ii);
  leg{ii}=nemo_t(kk).name;
end
pwd_old=pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ MOMENTS AND ASSOCIATED DATA FOR TRANSPORT %%%%%%%%%%%%%%%%%%%%
for jj=1:length(ind)
  %
  k=ind(jj);
  cd(nemo_t(k).path);
  filename=nemo_t(k).filename;
  for ii=1:size(nemo_t(k).kineticspecies,2)
    s0d=nemo_t(k).kineticspecies{ii};
    if strcmp(nemo_t(k).kineticspecies{ii},'generic') == 0
      disp(s0d)
      if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	s1=2;
	s2=1;
	nemo_t(k).irad_coord=2;
      else
	s1=1;
	s2=0;
	nemo_t(k).irad_coord=1;
      end
      %
      nemo_t(k).(s0d).time_1D    = double(hdf5read(filename,['/data/var1d/',s0d,'/f_av/time']));
      tmp = double(hdf5read(filename,['/data/var1d/',s0d,'/efluxw_rad/coord1']));
      nemo_t(k).(s0d).psi_prof1D = tmp(s1:end-s2);
      nemo_t(k).(s0d).s_prof1D 	 = sqrt(nemo_t(k).(s0d).psi_prof1D);
      nemo_t(k).(s0d).nprof_1D 	 = size(nemo_t(k).(s0d).v_perp2_av,1);
      nemo_t(k).(s0d).nsteps_1D	 = size(nemo_t(k).(s0d).v_perp2_av,2);
      %
      % Get d(rho/a)/ds on a \tilde(psi) grid
      %
      if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	for ll=1:nemo_t(k).(s0d).nprof_1D-(s1+s2-1)
	  %get rho/a
	  psi_index=find(nemo_t(k).(s0d).psi_prof < ...
			 nemo_t(k).(s0d).psi_prof1D(ll),1,'last');
	  if isempty(psi_index) 
	    psi_index=1;
	  end
	  if psi_index > size(nemo_t(k).(s0d).psi_prof,1)-1
	    psi_index=size(nemo_t(k).(s0d).psi_prof,1)-1;
	    w1=1;
	  else
	    w1=(nemo_t(k).(s0d).psi_prof1D(ll)-nemo_t(k).(s0d).psi_prof(psi_index))/(nemo_t(k).(s0d).psi_prof(psi_index+1)-nemo_t(k).(s0d).psi_prof(psi_index));
	  end
	  w0=1-w1;
	  psi_index;
	  %rho_int=rho/a (psi_prof1D(i))
	  rho_int=w0*nemo_t(k).(s0d).rho_prof(psi_index)+w1*nemo_t(k).(s0d).rho_prof(psi_index+1);
	  nemo_t(k).vol_prof1D(ll)=rho_int;
	  rho_index=find(nemo_t(k).(s0d).rho_prof< rho_int,1,'last');
	  if isempty(rho_index)
	    rho_index=1;
	  end
	  if rho_index > size(nemo_t(k).(s0d).rho_prof)-1
	    rho_index=size(nemo_t(k).(s0d).rho_prof)-1;
	  end
	  w1=(rho_int-nemo_t(k).(s0d).rho_prof(rho_index))/(nemo_t(k).(s0d).rho_prof(rho_index+1)-nemo_t(k).(s0d).rho_prof(rho_index));
	  w0=1-w1;
	  nemo_t(k).(s0d).drhods_int(ll)=w0*nemo_t(k).(s0d).drhods(rho_index)+w1*nemo_t(k).(s0d).drhods(rho_index+1);
	end

	if size(nemo_t(k).(s0d).s_prof,1) < size(nemo_t(k).(s0d).s_prof,2)
	  nemo_t(k).(s0d).s_prof=nemo_t(k).(s0d).s_prof';
	end

	%select radial coordinates for plot
	switch nemo_t(k).irad_coord
	 case 1
	  nemo_t(k).radlab='s';
	  nemo_t(k).(s0d).rad_prof1D=nemo_t(k).(s0d).s_prof1D;
	  nemo_t(k).(s0d).rad_prof=nemo_t(k).(s0d).s_prof;
	 case 2
	  nemo_t(k).radlab='sqrt(V(s)/V)';
	  nemo_t(k).(s0d).rad_prof1D=nemo_t(k).vol_prof1D;
	  nemo_t(k).(s0d).rad_prof=nemo_t(k).(s0d).rho_prof;
	end
	%
% $$$ elseif nemo_t(k).(s0d).nsel_profile == 5
% $$$ nemo_t(k).radlab='s';
% $$$ nemo_t(k).(s0d).rad_prof1D=nemo_t(k).(s0d).s_prof1D;
% $$$ ngmax=size(nemo_t(k).(s0d).t_pic,1);
% $$$ dpsi_0 = 1/ngmax;
% $$$ psi_grd = [0:dpsi_0:1-dpsi_0];
% $$$ nemo_t(k).(s0d).rad_prof=sqrt(psi_grd) ; % sc.coord1
% $$$ nemo_t(k).(s0d).s_prof=sqrt(psi_grd)';
      else %adhoc && nsel_profile
	nemo_t(k).radlab='s';
	nemo_t(k).(s0d).rad_prof1D=nemo_t(k).(s0d).s_prof1D;
	nemo_t(k).(s0d).rad_prof=nemo_t(k).(s0d).s_prof;
      end %adhoc && nsel_profile

      n=nemo_t(k).(s0d).nsteps_1D;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE IONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      nemo_t(k).(s0d).mytemp=1./3*nemo_t(k).(s0d).mass* ...
	  (nemo_t(k).(s0d).v_perp2_av./nemo_t(k).(s0d).f_av+nemo_t(k).(s0d).v_par2_av./nemo_t(k).(s0d).f_av-(nemo_t(k).(s0d).v_par_av./nemo_t(k).(s0d).f_av).^2);
      %
      % Smoothing (cubic splines) of temperature and temperature gradient 
      %
      disp(['Computing ', s0d, ' temperature gradient for ', nemo_t(k).name])
      for ll=2:n
	int_mytemp=csaps(double(nemo_t(k).(s0d).psi_prof1D),double(nemo_t(k).(s0d).mytemp(s1:end-s2,ll)),1.0-tol);
	%polynomial coefficient of temperature gradient
	dint_mytemp=fnder(int_mytemp);
	%dT/d\tilde(psi) 
	dmytemp=ppval(nemo_t(k).(s0d).psi_prof1D,dint_mytemp);
	%dT/ds, normalised on a higher resolution \tilde(psi) grid
	%Multiply by tau to get the correct normalization
	nemo_t(k).(s0d).dTds(:,ll)=dmytemp.*nemo_t(k).(s0d).s_prof1D*2*nemo_t(k).(s0d).tau;
      end
      % Get dT(rho/a)/d(rho/a)
      if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile ...
	    > 1
	for ll=2:n
	  nemo_t(k).(s0d).dTdrhooa(:,ll)=nemo_t(k).(s0d).dTds(:,ll)./transpose(nemo_t(k).(s0d).drhods_int);
	end
      end
      %
      if lplot
	%plot analytic and reconstructed temperature
	figure;
	hold on
	plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).t_pic,'b')
	plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).mytemp(s1:end-s2,1),'r')
	plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).mytemp(s1:end-s2,n),'g')
	axis ([0 1 0 max(nemo_t(k).(s0d).t_pic)+1]);
	xlabel(nemo_t(k).radlab)
	ylabel('T_i/T_i(s_0)')
	legend('analytic','particles, initial','particles, final')
	title(strcat(s0d, ' temperature profile, ',nemo_t(k).name))

	figure;
	hold on

	if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	  %plot initial analytic and reconstructed logarithmic temperature gradient

	  plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).gradt_pic,'b')
	  plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dTdrhooa(:,2)./nemo_t(k).(s0d).mytemp(s1:end-s2,2),'r')
	  plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dTdrhooa(:,n)./nemo_t(k).(s0d).mytemp(s1:end-s2,n),'k')
	  ylabel('a/L_{T,i,\rho/a}')
	else
	  % Normalise analytic temperature gradient
	  nemo_t(k).(s0d).gradt_pic=nemo_t(k).(s0d).gradt_pic.*nemo_t(k).(s0d).s_prof*2*nemo_t(k).(s0d).tau;
	  plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).gradt_pic,'b')
	  plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dTds(:,2)./nemo_t(k).(s0d).mytemp(s1:end-s2,2),'r')
	  plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dTds(:,n)./nemo_t(k).(s0d).mytemp(s1:end-s2,n),'k')
	  ylabel('a/L_{T,i,s}')
	end

	axis ([0 1 (-nemo_t(k).(s0d).kappat0)-1 max(nemo_t(k).(s0d).gradt_pic)+1]);
	xlabel(nemo_t(k).radlab)
	legend('analytic','particles, initial','particle, final')
	title(strcat(s0d, ' temperature gradient profile, ',nemo_t(k).name))
      end %lplot
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  disp('Smoothing density and gradient')
	  for ll=2:n
	    int_f_av=csaps(double(nemo_t(k).(s0d).psi_prof1D),double(nemo_t(k).(s0d).f_av(s1:end-s2,ll)),1.0-tol);
	    %polynomial coefficient of density gradient
	    dint_f_av=fnder(int_f_av);
	    %dn/d\tilde(psi) on a \tilde(psi) grid
	    df_av=ppval(nemo_t(k).(s0d).psi_prof1D,dint_f_av);
	    %dn/ds, normalised on a higher resolution \tilde(psi) grid
	    nemo_t(k).(s0d).dnds(:,ll)=df_av.*nemo_t(k).(s0d).s_prof1D*2;         
	  end
	  %Get dn/d(rho/a) on a \tilde(psi) grid
	  if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	    for ll=2:n
	      nemo_t(k).(s0d).dndrhooa(:,ll)=nemo_t(k).(s0d).dnds(:,ll)./transpose(nemo_t(k).(s0d).drhods_int);
	    end
	  end
	  %
	  if lplot
	    %plot analytic and reconstructed density
	    figure;
	    hold on
	    plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).n_pic,'b')
	    plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).f_av(s1:end-s2,2),'r')
	    plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).f_av(s1:end-s2,n),'g')
	    axis ([0 1 0 max(nemo_t(k).(s0d).n_pic)+1]);
	    xlabel(nemo_t(k).radlab)
	    ylabel('n_i/<n_i>')
	    legend('analytic','particles, initial','particles, final')
	    title(strcat('Ion density profile, ',nemo_t(k).name))
	    %plot analytic and reconstructed logarithmic density gradient
	    figure;
	    hold on
	    if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	      plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).gradn_pic,'b')
	      plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dndrhooa(:,2)./nemo_t(k).(s0d).f_av(s1:end-s2,2),'r')
	      plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dndrhooa(:,n)./nemo_t(k).(s0d).f_av(s1:end-s2,n),'k')
	      ylabel('a/L_{n,i,\rho/a}')
	    else
	      % Normalise analytic density gradient
	      nemo_t(k).(s0d).gradn_pic=nemo_t(k).(s0d).gradn_pic.*nemo_t(k).(s0d).s_prof*2;
	      plot(nemo_t(k).(s0d).rad_prof,nemo_t(k).(s0d).gradn_pic,'b')
	      plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dnds(:,2)./nemo_t(k).(s0d).f_av(s1:end-s2,1),'r')
	      plot(nemo_t(k).(s0d).rad_prof1D,nemo_t(k).(s0d).dnds(:,n)./nemo_t(k).(s0d).f_av(s1:end-s2,n),'k')
	      ylabel('a/L_{n,i,s}')
	    end
	    axis ([0 1 (-nemo_t(k).(s0d).kappan0)-1 max(nemo_t(k).(s0d).gradn_pic)+1]);
	    xlabel(nemo_t(k).radlab)
	    legend('analytic','particles, initial','particles, final')
	    title(strcat('Ion density gradient profile, ',nemo_t(k).name))
	  end

	  if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	    nemo_t(k).(s0d).RoLT=-nemo_t(k).(s0d).dTdrhooa./nemo_t(k).(s0d).mytemp(s1:end-s2,:)*nemo_t(k).generic.aspect_ratio;
	  else
	    nemo_t(k).(s0d).RoLT=-nemo_t(k).(s0d).dTds./nemo_t(k).(s0d).mytemp(s1:end-s2,:)*nemo_t(k).generic.aspect_ratio;
	  end
	  %
	  if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1
	    nemo_t(k).(s0d).RoLn=-nemo_t(k).(s0d).dndrhooa./nemo_t(k).(s0d).f_av(s1:end-s2,:)*nemo_t(k).generic.aspect_ratio;
	  else
	    nemo_t(k).(s0d).RoLn=-nemo_t(k).(s0d).dnds./nemo_t(k).(s0d).f_av(s1:end-s2,:)*nemo_t(k).generic.aspect_ratio;
	  end

	  %%%%%%%%%%%%%%%%%%%%%%%%%%% ETAI IONS%%%%%%%%%%%%%%%%%%%%%
	  nemo_t(k).(s0d).etai=nemo_t(k).(s0d).RoLT./nemo_t(k).(s0d).RoLn;
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  %initialize chi_perp
	  nemo_t(k).(s0d).chi_perp=0*nemo_t(k).(s0d).efluxw_rad;
	  if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1 & nemo_t(k).norm_Dimits==1
	    nemo_t(k).(s0d).chi_perp=-nemo_t(k).(s0d).efluxw_rad(s1:end-s2,:)./(nemo_t(k).(s0d).f_av(s1:end-s2,:).*nemo_t(k).(s0d).dTdrhooa(:,:))*(nemo_t(k).generic.lx/2)^2/nemo_t(k).(s0d).kappan0;
	  else
	    nemo_t(k).(s0d).chi_perp=-nemo_t(k).(s0d).efluxw_rad(s1:end-s2,:)./(nemo_t(k).(s0d).f_av(s1:end-s2,:).*nemo_t(k).(s0d).dTds(:,:))*(nemo_t(k).generic.lx/2)^2;
	  end

	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  %initialize D_perp
	  nemo_t(k).(s0d).D_perp=0*nemo_t(k).(s0d).pfluxw_rad(s1:end-s2);
	  if nemo_t(k).generic.nsel_equil == 2 & nemo_t(k).(s0d).nsel_profile > 1 & nemo_t(k).norm_Dimits==1
	    nemo_t(k).(s0d).D_perp = -nemo_t(k).(s0d).pfluxw_rad(s1:end-s2,:)./...
		nemo_t(k).(s0d).dndrhooa(s1:end-2,:)*(nemo_t(k).generic.lx/2)^2/nemo_t(k).(s0d).kappan0;
	  else
	    nemo_t(k).(s0d).D_perp =-nemo_t(k).(s0d).pfluxw_rad(s1:end-s2,:)./nemo_t(k).(s0d).dnds(:,:)*(nemo_t(k).generic.lx/2)^2;
	  end

    end
  end %loop over nemo_
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %species loop
nemo_t=nemo_t;
cd(pwd_old);

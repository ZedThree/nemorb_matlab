function omegar=nemorb_linear_omega(sim,s0, s0d, lplot,ind)
%######################################
%omegar=real_freq_hdf5(sim,s0, s0d, lplot,ind)
%######################################
%---------------
%
%Description
%
%For linear runs, compute real frequency of toroidal mode
%
%---------------
%
%Dummy arguments
%
%sim   = array containing the simulations data, created with loader.m
%s0    = value of magnetic surface to extract the real frequency
%s0d   = species name
%lplot : true = plot results
%ind   = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%omegar = real frequency(1) + left(2) and right(3) values in \Omega_i^{-1} units
%---------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REAL FREQUENCY (LINEAR SIMULATIONS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  global timelab

  if exist('ind')==0
    ind=1;
  end
  colorp={'blue','red','green','black','cyan','yellow'};

  pwd_old=pwd;

  for ii=1:length(ind) %loop on simulations
    k=ind(ii);
    %go to simulation directory
    cd(sim(k).path);

    %check that tmin_rf and tmax_rf are correct
    if sim(k).tmax_rf > sim(k).time(end)
      sim(k).tmax_rf=sim(k).time(end);
      disp(sprintf('tmax_rf modified because it was too big, new tmax_rf = %15.5f', sim(k).tmax_rf))
    end
    if sim(k).tmin_rf > sim(k).time(end)
      sim(k).tmin_rf=sim(k).time(1);
      disp(sprintf('tmin_rf modified because it was too big, new tmin_rf = %15.5f', sim(k).tmin_rf))
    end

    imin_rf = find(sim(k).(s0d).time_1D >= sim(k).tmin_rf,1,'first');
    imax_rf = find(sim(k).(s0d).time_1D <= sim(k).tmax_rf,1,'last');

    %number of discrete points for the fft
    npoints_fft=imax_rf-imin_rf+1;

    %modify imin_rf to have an even number of points
    if mod(npoints_fft,2)==1
      imax_rf=imax_rf-1;
      disp('tmax_rf modified to have even number of points for fft')
      npoints_fft=imax_rf-imin_rf+1;
    end

    %read data
    s_freq=hdf5read(sim(k).filename, '/data/var1d/generic/philoc/coord1');
    phi_freq=hdf5read(sim(k).filename, '/data/var1d/generic/philoc/data');
    is_freq=find(s_freq <= s0,1,'last');
    if s_freq(is_freq) ~= s0
      disp(sprintf('Data does not exist for s = %0.5g',s0))
      sreal=s_freq(is_freq);
      disp(sprintf('extract real frequency at s = %0.5g',sreal))
    end


    phi_fit=transpose(squeeze(double(log(abs(phi_freq(is_freq,imin_rf:imax_rf))))));
    time_fit=double(sim(k).(s0d).time_1D(imin_rf:imax_rf));
    %extract a growth rate
    [coef,fit]=polyfit(time_fit,phi_fit,1);
    gamma_loc(ii)=coef(1);
    disp(sprintf('Growth rate = %0.5e',gamma_loc(ii)))
    if(lplot)
    figure
    plot(sim(k).(s0d).time_1D,double(log(abs(phi_freq))))
    xlabel('s')
    ylabel('time')
    end
    %modify potential;
    phi0m=transpose(phi_freq(is_freq,:)).*exp(-gamma_loc(ii)*sim(k).(s0d).time_1D);
    phi0m=double(phi0m);

    delf=2*pi/sim(k).generic.dt;
    %range of frequency
    harmonics=[-npoints_fft/2:npoints_fft/2-1];
    freq=delf*harmonics;

    %DFT, shifted to the center of the spectrum
    phifft=fftshift(fft(phi0m(imin_rf:imax_rf)));

    %get maximum component and its location
    [phimax_rf,imax_rf_freq]=max(abs(real(phifft)));
    %get maximum frequency
    fmax=freq(imax_rf_freq)/npoints_fft;

    %Linear interpolations to have a better approximation in
    %case there is a lack of resolution
    fmax_l=freq(imax_rf_freq-1)/npoints_fft;
    fmax_r=freq(imax_rf_freq+1)/npoints_fft;

    phimax_rf_l=abs(real(phifft(imax_rf_freq-1)));
    phimax_rf_r=abs(real(phifft(imax_rf_freq+1)));
    phimax_rf_sum_l=phimax_rf_l+phimax_rf;
    phimax_rf_sum_r=phimax_rf_r+phimax_rf;

    wl0=phimax_rf_l/phimax_rf_sum_l;
    wl1=1-wl0;
    wr0=phimax_rf/phimax_rf_sum_r;
    wr1=1-wr0;

    fmax_l_int=wl0*fmax_l+wl1*fmax;
    fmax_r_int=wr0*fmax+wr1*fmax_r;

    disp(sprintf('Frequency of the mode = %0.5e',    fmax))
    disp(sprintf('Frequency of the mode, left = %0.5e',  fmax_l_int))
    disp(sprintf('Frequency of the mode, right = %0.5e', fmax_r_int))

    %store 3 frequencies
    omegar(ii,1)=fmax;
    omegar(ii,2)=fmax_l_int;
    omegar(ii,3)=fmax_r_int;

    %Try to fit a cosinus function 
    phifft_fit=phifft;
    %select only the main 2 components
    for j=1:length(phifft)
      if  xor(j~=imax_rf_freq,j~=npoints_fft+2-imax_rf_freq)==0
	phifft_fit(j)=complex(0,0);
      end   
    end

    %inverse Fourier transform (imaginary part is very small, so we take the real part)
    phi_fit=real(ifft(ifftshift(phifft_fit))); %phi_fit=cosinus function

    %read phi(s_0, chi,t)
    phi0_chi=hdf5read(sim(k).filename,'/data/var1d/generic/phi0_chi/data');
    chi_1D=hdf5read(sim(k).filename,'/data/var1d/generic/phi0_chi/coord1');
    %normalise phi0_chi
    phi0_chi_norm=phi0_chi;
    for j=1:length(sim(k).(s0d).time_1D)
      phi0_chi_norm(:,j)=phi0_chi(:,j)/max(abs(phi0_chi(:,j)),[],1);
    end

    if(lplot)
      %plot abs(phi0) on a log scale
      figure;
      hold on
      plot(sim(k).(s0d).time_1D,log(abs(phi_freq(is_freq,:))));
      plot(sim(k).(s0d).time_1D(imin_rf:imax_rf),coef(2)+coef(1)*sim(k).(s0d).time_1D(imin_rf:imax_rf),'r-')
      xlabel(timelab{1})
      ylabel('\phi(s_0,0,0) [q_i/T_e(s_0)]')
      title(sim(k).name)
      %plot modified potential
      figure;
      plot(sim(k).(s0d).time_1D,phi0m)
      hold on
      xlabel(timelab{1})
      ylabel('\phi(s_0,0,0)e^{-\gamma t}')
      plot(sim(k).(s0d).time_1D(imin_rf:imax_rf),phi_fit,'r')
      title(sim(k).name)
      legend('ORB5','fit')

      %plot frequency spectrum around maximum
      i1=max(1,imax_rf_freq-sim(k).wfs);
      i2=min(npoints_fft,imax_rf_freq+sim(k).wfs);
      figure;
      plot(freq(i1:i2),real(phifft(i1:i2)),'bx-');
      grid on;
      legend('real part')
      title(strcat(sim(k).name,': estimated frequency \omega=',num2str(abs(fmax))))

      %plot phi(s0,chi,0,t) to have the sign of the frequency
      figure;
      pcolor(sim(k).(s0d).time_1D,double(chi_1D),double(phi0_chi_norm))
      shading interp;
      colorbar
      xlabel(timelab{1})
      if sim(k).generic.nsel_coord==1
	ylabel('\theta')
	title(strcat('\phi(s_0,\theta,0,t), ',sim(k).name))
      else
	ylabel('\theta_*')
	title(strcat('\phi(s_0,\theta_*,0,t), ',sim(k).name))
      end
    end %lplot
  end %loop on simulations

  cd(pwd_old)

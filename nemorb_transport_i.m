function nemorb_transport_i(sim, ind)
%################################
%nemorb_transport_i(sim, ind)
%################################
%---------------
%
%Description
%
%Do all the transport and suchlike
%
%---------------
%
%Dummy arguments
%
%sim = array containing the simulations data, created with loader.m
%ind = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%---------------

global timelab
if exist('ind') == 0
ind=1;
end

pwd_old=pwd;
for i=1:length(ind) %Loop on simulations
k=ind(i);
cd(sim(k).path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE IONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mytemp=1./3*(struct_out(i).v_perp2_av ./ struct_out(i).f_av + ...
	     struct_out(i).v_par2_av ./ struct_out(i).f_av - ...
	     (struct_out(i).v_par_av ./ struct_out(i).f_av).^2);
%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, ions')
for i=2:nsteps_1D
int_mytemp=csaps(double(struct_out(i).psi_prof1D),double(mytemp(:,i)),1-tol);
%coefficient of temperature gradient
dint_mytemp=fnder(int_mytemp);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i=ppval(struct_out(i).psi_prof1D_hr,dint_mytemp);

% PH: fixing 1D_hr (nb: hr = higher resolution?)

%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i(:,i)=ppval(struct_out(i).psi_prof1D_hr,int_mytemp);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds(:,i)=dmytemp_i.*s_prof1D_hr*2;
end

%
% Get dT(rho/a)/d(rho/a) on a higher resolution \tilde(psi) grid
%

if adhoc ==1 & nsel_profile > 1
for i=2:nsteps_1D
dTdrhooa(:,i)=dTds(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed temperature
figure;
plot(s_prof,ti,'b',s_prof1D,mytemp(:,1),'r',s_prof1D_hr,mytemp_i(:,1),'k',s_prof1D,mytemp(:,nsteps_1D),'g',s_prof1D_hr,mytemp_i(:,nsteps_1D),'c')
xlabel('s')
ylabel('T_i/T_e(s_0)')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Initial ion temperature profile')

figure;
if adhoc == 1 & nsel_profile >1
%plot initial analytic and reconstructed logarithmic temperature gradient
plot(s_prof,gradt,'b',s_prof1D_hr,dTdrhooa(:,1)./mytemp_i(:,1),'r',s_prof1D_hr,dTdrhooa(:,nsteps_1D)./mytemp_i(:,nsteps_1D),'k')
ylabel('a/L_{T,i,\rho/a}')
else
plot(s_prof,gradt,'b',s_prof1D_hr,dTds(:,1)./mytemp_i(:,1),'r',s_prof1D_hr,dTds(:,nsteps_1D)./mytemp_i(:,nsteps_1D),'k')
ylabel('a/L_{T,i,s}')
end
xlabel('s')
legend('analytic','particles, initial','particle, final')
title('Ion temperature gradient profile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Smoothing density and gradient')
for i=2:nsteps_1D
int_f_av=csaps(double(psi_prof1D),double(f_av(:,i)),1-tol);
f_av_i(:,i)=ppval(psi_prof1D_hr,int_f_av);
%coefficient of density gradient
dint_f_av=fnder(int_f_av);
%dn/d\tilde(psi) on a higher resolution \tilde(psi) grid
df_av_i=ppval(psi_prof1D_hr,dint_f_av);
%dn/ds, normalised on a higher resolution \tilde(psi) grid
dnds(:,i)=df_av_i.*s_prof1D_hr*2;         
end
%
% Get dn/d(rho/a) on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=2:nsteps_1D
dndrhooa(:,i)=dnds(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed density
figure;
plot(s_prof,ni,'b',s_prof1D,f_av(:,1),'r',s_prof1D_hr,f_av_i(:,1),'k',s_prof1D,f_av(:,nsteps_1D),'g',s_prof1D_hr,f_av_i(:,nsteps_1D),'c')
xlabel('s')
ylabel('n_i/<n_i>')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Ion density profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradn,'b',s_prof1D_hr,dndrhooa(:,1)./f_av_i(:,1),'r',s_prof1D_hr,dndrhooa(:,nsteps_1D)./f_av_i(:,nsteps_1D),'k')
ylabel('a/L_{n,i,\rho/a}')
else
plot(s_prof,gradn,'b',s_prof1D_hr,dnds(:,1)./f_av_i(:,1),'r',s_prof1D_hr, dnds(:,nsteps_1D)./f_av_i(:,nsteps_1D),'k')
ylabel('a/L_{n,i,s}')
end
xlabel('s')
legend('analytic','particles, initial','particles, final')
title('Ion density gradient profile')
if adhoc == 1 & nsel_profile >1
RoLT=-dTdrhooa./mytemp_i*aspect_ratio;
else
RoLT=-dTds./mytemp_i*aspect_ratio;
end


%
%Plot evolution of R/Ln
%
if adhoc == 1 & nsel_profile >1
RoLn=-dndrhooa./f_av_i*aspect_ratio;
else
RoLn=-dnds./f_av_i*aspect_ratio;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAI IONS%%%%%%%%%%%%%%%%%%%%%
etai=RoLT./RoLn;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM IONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing ion parallel momentum')
for i=2:nsteps_1D
int_v_par_av=csaps(double(psi_prof1D),double(v_par_av(:,i)),1-tol);
%v_par_av on a higher resolution \tilde(psi) grid
v_par_av_hr(:,i)=ppval(psi_prof1D_hr,int_v_par_av);
end

v_par_av_hr=v_par_av_hr./f_av_i;
v_par_av_hr=double(v_par_av_hr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing ion heatflux')
for i=2:nsteps_1D
int_efluxw_rad=csaps(double(psi_prof1D),double(efluxw_rad(:,i)),1-tol);
efluxw_rad_i(:,i)=ppval(psi_prof1D_hr,int_efluxw_rad);
end

%initialize chi_perp
chi_perp=0*efluxw_rad_i;
if adhoc == 1 & nsel_profile >1 & norm_Dimits==1
chi_perp(imin_1D:imax_1D,:)=-efluxw_rad_i(imin_1D:imax_1D,:)./(f_av_i(imin_1D:imax_1D,:).*dTdrhooa(imin_1D:imax_1D,:))*(lx/2)^2/kappan0;
else
chi_perp(imin_1D:imax_1D,:)=-efluxw_rad_i(imin_1D:imax_1D,:)./(f_av_i(imin_1D:imax_1D,:).*dTds(imin_1D:imax_1D,:))*(lx/2)^2;
end

%Plot chi vs R/LT
ilow=find(rad_prof1D_hr < radlow,1,'last');
iup =find(rad_prof1D_hr < radup,1,'last');
text=sprintf('%s', 'Profiles are averaged between ',radlab,'=');
disp(text);
rad_prof1D_hr(ilow)
rad_prof1D_hr(iup)

RoLT_av=squeeze(mean(RoLT(ilow:iup,:)));
chi_perp_av=squeeze(mean(chi_perp(ilow:iup,:)));
RoLT_av_loc(1)=RoLT_av(1);
chi_perp_av_loc(1)=chi_perp_av(1);

npoints=floor(nsteps_1D/tlength)-1;
for i=1:npoints
RoLT_av_loc(i+1)=RoLT_av(1+i*tlength);
chi_perp_av_loc(i+1)=chi_perp_av(1+i*tlength);
end

RoLT_graph=[0.8*min(RoLT_av_loc):0.01:max(4.1,1.2*max(RoLT_av_loc))];
if norm_Dimits==1
chi_perp_Dimits=15.4*(1-6./RoLT_graph);
else
chi_perp_Dimits=15.4*kappan0*(1-6./RoLT_graph);
end

figure;
hold on;
plot(RoLT_av_loc,chi_perp_av_loc,'kx-');
plot(RoLT_graph,chi_perp_Dimits,'b');
axis([RoLT_graph(1) RoLT_graph(end) 0 1.2*max(abs(chi_perp_av_loc))])
xlabel('R_0/L_{T,i}')
ylabel('\chi_i/\chi_{GB}')

figure;
hold on
plot(time_1D, RoLT_av,'k')
xlabel(timelab{1})
ylabel('R/L_{T,i}')

figure;
hold on
plot(time_1D, chi_perp_av,'k')
xlabel(timelab{1})
ylabel('\chi/\chi_{GB}')


end %simulation loop
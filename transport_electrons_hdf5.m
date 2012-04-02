%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE ELECTRONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlelec > 0

%get total contribution
f_av_el=f_av_p_el+f_av_tr_el;
v_par_av_el=v_par_av_p_el+v_par_av_tr_el;
v_par2_av_el=v_par2_av_p_el+v_par2_av_tr_el;
v_perp2_av_el=v_perp2_av_p_el+v_perp2_av_tr_el;
efluxw_rad_el=efluxw_rad_p_el+efluxw_rad_tr_el;
efluxf0_rad_el=efluxf0_rad_p_el+efluxf0_rad_tr_el;

mytemp_el    =melec/mion/3*(v_perp2_av_el./f_av_el+v_par2_av_el./f_av_el-(v_par_av_el./f_av_el).^2);
mytemp_el_nop =melec/mion/3*(v_perp2_av_tr_el./f_av_tr_el+v_par2_av_tr_el./f_av_tr_el-(v_par_av_tr_el./f_av_tr_el).^2);

%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, electrons')
for i=2:nsteps_1D
int_mytemp_el=csaps(double(psi_prof1D),double(mytemp_el(:,i)),1-tol);
int_mytemp_el_nop=csaps(double(psi_prof1D),double(mytemp_el_nop(:,i)),1-tol);
%coefficient of temperature gradient
dint_mytemp_el=fnder(int_mytemp_el);
dint_mytemp_el_nop=fnder(int_mytemp_el_nop);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i_el=ppval(psi_prof1D_hr,dint_mytemp_el);
dmytemp_i_el_nop=ppval(psi_prof1D_hr,dint_mytemp_el_nop);
%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i_el(:,i)=ppval(psi_prof1D_hr,int_mytemp_el);
mytemp_i_el_nop(:,i)=ppval(psi_prof1D_hr,int_mytemp_el_nop);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds_el(:,i)=dmytemp_i_el.*s_prof1D_hr*2;
dTds_el_nop(:,i)=dmytemp_i_el_nop.*s_prof1D_hr*2;
end
%
% Get dT/d(rho/a) on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=2:nsteps_1D
dTdrhooa_el(:,i)=dTds_el(:,i)./transpose(drhods_int);
dTdrhooa_el_nop(:,i)=dTds_el_nop(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed temperature
figure;
plot(s_prof,te,'b',s_prof1D,mytemp_el(:,2),'r',s_prof1D_hr,mytemp_i_el(:,2),'k',s_prof1D,mytemp_el(:,nsteps_1D),'g',s_prof1D_hr,mytemp_i_el(:,nsteps_1D),'c')
xlabel('s')
ylabel('T_e/T_e(s_0)')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Electron temperature profile')

%plot analytic, trapped and passing temperature
figure
plot(s_prof,te,'b', s_prof1D_hr,mytemp_i_el(:,nsteps_1D), 'r', s_prof1D_hr, mytemp_i_el_nop(:,nsteps_1D),'k')
xlabel('s')
ylabel('T_e/T_e(s_0)')
legend('analytic','with detrapped, final','without detrapped, final')
title('Passing and trapped electrons temperature profile')

%plot initial analytic and reconstructed logarithmic temperature gradient
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradte,'b',s_prof1D_hr,dTdrhooa_el(:,2)./mytemp_i_el(:,2),'r',s_prof1D_hr,dTdrhooa_el(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'k')
ylabel('a/L_{T,e,\rho/a}')
else
plot(s_prof,gradte,'b',s_prof1D_hr,dTds_el(:,2)./mytemp_i_el(:,2),'r',s_prof1D_hr,dTds_el(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'k')
ylabel('a/L_{T,e,s}')
end
xlabel('s')
legend('analytic','particles, initial','particle, final')
title('Electron temperature gradient profile')

%plot initial analytic and reconstructed logarithmic temperature gradient, with and without detrapped
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradte,'b',s_prof1D_hr,dTdrhooa_el(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'r',s_prof1D_hr,dTdrhooa_el_nop(:,nsteps_1D)./mytemp_i_el_nop(:,nsteps_1D),'k')
ylabel('a/L_{T,e,\rho/a}')
else
plot(s_prof,gradte,'b',s_prof1D_hr,dTds_el(:,nsteps_1D)./mytemp_i_el(:,2),'r',s_prof1D_hr,dTds_el_nop(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'k',s_prof1D_hr)
ylabel('a/L_{T,e,s}')
end
xlabel('s')
legend('analytic','total, final','trapped, final')
title('Electron temperature gradient profile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlelec==2
for i=2:nsteps_1D
f_av_el_tot(:,i)=f_av_el(:,i)+ne0p_1D;
f_av_el_tot_nop(:,i)=f_av_tr_el(:,i)+ne0p_1D;
end
else
f_av_el_tot=f_av_el;
f_av_el_tot_nop=f_av_el_tot;
end

disp('Smoothing electron density and gradient')
for i=2:nsteps_1D
int_f_av_el=csaps(double(psi_prof1D),double(f_av_el_tot(:,i)),1-tol);
int_f_av_el_nop=csaps(double(psi_prof1D),double(f_av_el_tot_nop(:,i)),1-tol);
f_av_i_el_tot(:,i)=ppval(psi_prof1D_hr,int_f_av_el);
f_av_i_el_tot_nop(:,i)=ppval(psi_prof1D_hr,int_f_av_el_nop);
%coefficient of density gradient
dint_f_av_el=fnder(int_f_av_el);
dint_f_av_el_nop=fnder(int_f_av_el_nop);
%dn/d\tilde(psi) on a higher resolution \tilde(psi) grid
df_av_i_el=ppval(psi_prof1D_hr,dint_f_av_el);
df_av_i_el_nop=ppval(psi_prof1D_hr,dint_f_av_el_nop);
%dn/ds, normalised on a higher resolution \tilde(psi) grid
dnds_el(:,i)=df_av_i_el.*s_prof1D_hr*2;
dnds_el_nop(:,i)=df_av_i_el_nop.*s_prof1D_hr*2;
end
%
% Get dn/d(rho/a) on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=2:nsteps_1D
dndrhooa_el(:,i)=dnds_el(:,i)./transpose(drhods_int);
dndrhooa_el_nop(:,i)=dnds_el_nop(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed density
figure;
plot(s_prof,ne,'b',s_prof1D,f_av_el_tot(:,2),'r',s_prof1D_hr,f_av_i_el_tot(:,2),'k',s_prof1D,f_av_el_tot(:,nsteps_1D),'g',s_prof1D_hr,f_av_i_el_tot(:,nsteps_1D),'c')
xlabel('s')
ylabel('n_e/<n_e>')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Electron density profile')

figure;
plot(s_prof,ne,'b',s_prof1D_hr,f_av_i_el_tot(:,nsteps_1D),'r',s_prof1D_hr,f_av_i_el_tot_nop(:,nsteps_1D),'k')
xlabel('s')
ylabel('n_e/<n_e>')
legend('analytic','final','no detrapped, final')
title('Electron density profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradne,'b',s_prof1D_hr,dndrhooa_el(:,2)./f_av_i_el_tot(:,2),'r',s_prof1D_hr,dndrhooa_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'k')
ylabel('a/L_{n,e,\rho/a}')
else
plot(s_prof,gradne,'b',s_prof1D_hr,dnds_el(:,2)./f_av_i_el_tot(:,2),'r',s_prof1D_hr, dnds_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'k')
ylabel('a/L_{n,e,s}')
end
xlabel('s')
legend('analytic','particles, initial','particles, final')
title('Electron density gradient profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradne,'b',s_prof1D_hr,dndrhooa_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'r',s_prof1D_hr,dndrhooa_el_nop(:,nsteps_1D)./f_av_i_el_tot_nop(:,nsteps_1D),'k')
ylabel('a/L_{n,e,\rho/a}')
else
plot(s_prof,gradne,'b',s_prof1D_hr,dnds_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'r',s_prof1D_hr,dnds_el_nop(:,nsteps_1D)./f_av_i_el_tot_nop(:,nsteps_1D),'k')
ylabel('a/L_{n,e,s}')
end
xlabel('s')
legend('analytic','final','no detrapped, final')
title('Electron density gradient profile')

if adhoc == 1 & nsel_profile >1
RoLT_el=-dTdrhooa_el./mytemp_i_el*aspect_ratio;
RoLT_el_nop=-dTdrhooa_el_nop./mytemp_i_el_nop*aspect_ratio;
else
RoLT_el=-dTds_el./mytemp_i_el*aspect_ratio;
RoLT_el_nop=-dTds_el_nop./mytemp_i_el_nop*aspect_ratio;
end

%
%Plot evolution of R/Ln
%
if adhoc == 1 & nsel_profile >1
RoLn_el=-dndrhooa_el./f_av_i_el_tot*aspect_ratio;
RoLn_el_nop=-dndrhooa_el_nop./f_av_i_el_tot_nop*aspect_ratio;
else
RoLn_el=-dnds_el./f_av_i_el_tot*aspect_ratio;
RoLn_el_nop=-dnds_el_nop./f_av_i_el_tot_nop*aspect_ratio;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM ELECTRONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing electron parallel momentum')
for i=2:nsteps_1D
int_v_par_av_el=csaps(double(psi_prof1D),double(v_par_av_el(:,i)),1-tol);
int_v_par_av_el_nop=csaps(double(psi_prof1D),double(v_par_av_tr_el(:,i)),1-tol);
%v_par_av on a higher resolution \tilde(psi) grid
v_par_av_hr_el(:,i)=ppval(psi_prof1D_hr,int_v_par_av_el);
v_par_av_hr_el_nop(:,i)=ppval(psi_prof1D_hr,int_v_par_av_el_nop);
end

v_par_av_hr_el=v_par_av_hr_el./f_av_i_el_tot;
v_par_av_hr_el_nop=v_par_av_hr_el_nop./f_av_i_el_tot_nop;

%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAE ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etae=RoLT_el./RoLn_el;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing electron heatflux')
for i=2:nsteps_1D
int_efluxw_rad_el=csaps(double(psi_prof1D),double(efluxw_rad_el(:,i)),1-tol);
int_efluxw_rad_el_nop=csaps(double(psi_prof1D),double(efluxw_rad_tr_el(:,i)),1-tol);
efluxw_rad_i_el(:,i)=ppval(psi_prof1D_hr,int_efluxw_rad_el);
efluxw_rad_i_el_nop(:,i)=ppval(psi_prof1D_hr,int_efluxw_rad_el_nop);
end
%
%Plot evolution of Q
%

%initialize chi_perp
chi_perp_el=0*efluxw_rad_i_el;
chi_perp_el_nop=0*efluxw_rad_i_el;

chi_perp_el(imin_1D:imax_1D,:)=-efluxw_rad_i_el(imin_1D:imax_1D,:)./(f_av_i_el_tot(imin_1D:imax_1D,:).*dTdrhooa_el(imin_1D:imax_1D,:))*(lx/2)^2/kappan0;
chi_perp_el_nop(imin_1D:imax_1D,:)=-efluxw_rad_i_el_nop(imin_1D:imax_1D,:)./(f_av_i_el_tot_nop(imin_1D:imax_1D,:).*dTdrhooa_el_nop(imin_1D:imax_1D,:))*(lx/2)^2/kappan0;

if adhoc == 1 & nsel_profile >1 & norm_Dimits==1
else
chi_perp_el=chi_perp_el*kappan0;
chi_perp_el_nop=chi_perp_el_nop*kappan0;
end

%Plot chi vs R/LT
ilow=find(rad_prof1D_hr < radlow,1,'last');
iup =find(rad_prof1D_hr < radup,1,'last');
text=sprintf('%s', 'Profiles are averaged between ',radlab,'=');
disp(text);
rad_prof1D_hr(ilow)
rad_prof1D_hr(iup)

RoLT_av_el=squeeze(mean(RoLT_el(ilow:iup,:)));
RoLT_av_el_nop=squeeze(mean(RoLT_el_nop(ilow:iup,:)));
chi_perp_av_el=squeeze(mean(chi_perp_el(ilow:iup,:)));
chi_perp_av_el_nop=squeeze(mean(chi_perp_el_nop(ilow:iup,:)));

RoLT_av_loc_el(1)=RoLT_av_el(1);
RoLT_av_loc_el_nop(1)=RoLT_av_el_nop(1);
chi_perp_av_loc_el_nop(1)=chi_perp_av_el_nop(1);

npoints=floor(nsteps_1D/tlength)-1;
for i=1:npoints
RoLT_av_loc_el(i+1)=RoLT_av_el(1+i*tlength);
RoLT_av_loc_el_nop(i+1)=RoLT_av_el_nop(1+i*tlength);
chi_perp_av_loc_el(i+1)=chi_perp_av_el(1+i*tlength);
chi_perp_av_loc_el_nop(i+1)=chi_perp_av_el_nop(1+i*tlength);
end

RoLT_graph_el=[0.8*min(RoLT_av_loc_el):0.01:max(4.1,1.2*max(RoLT_av_loc_el))];
if norm_Dimits==1
chi_perp_Dimits_el=15.4*(1-6./RoLT_graph_el);
else
chi_perp_Dimits_el=15.4*kappan0*(1-6./RoLT_graph_el);
end

figure;
hold on;
plot(RoLT_av_loc_el,chi_perp_av_loc_el,'rx-');
plot(RoLT_av_loc_el_nop,chi_perp_av_loc_el_nop,'gx-');
if exist('chi_perp_av_loc') == 1
plot(RoLT_av_loc,chi_perp_av_loc,'bx-');
end
plot(RoLT_graph_el,chi_perp_Dimits_el,'k');
if exist('chi_perp_Dimits') == 1
plot(RoLT_graph,chi_perp_Dimits,'b');
end
if exist('chi_perp_av_loc') == 1
legend('electrons with detrapped','electrons without detrapped', 'ions')
else
legend('electrons with detrapped','electrons without detrapped')
end
xlabel('R_0/L_{T}')
ylabel('\chi/\chi_{GB}')

if exist('chi_perp_av_loc') == 1
chimaxi=max(chi_perp_av_loc);
RoLTmini=min(RoLT_av_loc);
RoLTmaxi=max(RoLT_av_loc);
else
chimaxi=-10000;
RoLTmini=10000;
RoLTmaxi=-10000;
end

chimaxe=max(chi_perp_av_loc_el);
chimax=max(chimaxi, chimaxe);
RoLTmine=min(RoLT_av_loc_el);
RoLTmin=min(RoLTmini, RoLTmine);
RoLTmaxe=max(RoLT_av_loc_el);
RoLTmax=max(RoLTmaxi, RoLTmaxe);

axis([0.8*RoLTmin 1.2*RoLTmax 0 1.2*chimax])

figure;
hold on
plot(time_1D, RoLT_av,'r',time_1D, RoLT_av_el,'b')
xlabel(timelab{1})
ylabel('R/L_{T,i}')
legend('ions', 'electrons')

figure;
hold on
plot(time_1D, chi_perp_av,'r',time_1D, chi_perp_av_el,'b')
xlabel(timelab{1})
ylabel('\chi/\chi_{GB}')
legend('ions','electrons')

end

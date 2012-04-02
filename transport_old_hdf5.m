%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE IONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mytemp=1./3*(v_perp2_av./f_av+v_par2_av./f_av-(v_par_av./f_av).^2);
psi_prof1D=coord1D;
psifmin=sfmin*sfmin;
psifmax=sfmax*sfmax;
psi_prof1D_hr=[psifmin:(psifmax-psifmin)/(nhr-1):psifmax];
s_prof1D=sqrt(psi_prof1D);
s_prof1D_hr=sqrt(psi_prof1D_hr);
s_prof1D_hr=double(s_prof1D_hr);


nprof_1D=size(v_perp2_av,1);
nsteps_1D=size(v_perp2_av,2);

%
% Get d(rho/a)/ds on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=1:nhr
%get rho/a
psi_index=find(psi_prof<psi_prof1D_hr(i),1,'last');
if isempty(psi_index) 
psi_index=1;
end
w1=(psi_prof1D_hr(i)-psi_prof(psi_index))/(psi_prof(psi_index+1)-psi_prof(psi_index));
w0=1-w1;
%rho_int=rho/a (psi_prof1D_hr(i))
rho_int=w0*rho_prof(psi_index)+w1*rho_prof(psi_index+1);
vol_prof1D_hr(i)=rho_int;
rho_index=find(rho_prof< rho_int,1,'last');
if isempty(rho_index)
rho_index=1;
end
w1=(rho_int-rho_prof(rho_index))/(rho_prof(rho_index+1)-rho_prof(rho_index));
w0=1-w1;
drhods_int(i)=w0*drhods(rho_index)+w1*drhods(rho_index+1);
end
end

%select radial coordinates for plot
switch irad_coord
case 1
rad_prof1D_hr=s_prof1D_hr;
radlab='s';
case 2
rad_prof1D_hr=vol_prof1D_hr;
radlab='sqrt(V(s)/V)';
end

imin_1D=find(rad_prof1D_hr < radmin_1D,1,'last');
imax_1D=find(rad_prof1D_hr < radmax_1D,1,'last');

%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, ions')
for i=1:nsteps_1D
int_mytemp=csaps(double(psi_prof1D),double(mytemp(:,i)),1-tol);
%coefficient of temperature gradient
dint_mytemp=fnder(int_mytemp);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i=ppval(psi_prof1D_hr,dint_mytemp);
%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i(:,i)=ppval(psi_prof1D_hr,int_mytemp);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds(:,i)=dmytemp_i.*s_prof1D_hr*2;
end

%
% Get dT(rho/a)/d(rho/a) on a higher resolution \tilde(psi) grid
%

if adhoc ==1 & nsel_profile > 1
for i=1:nsteps_1D
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
for i=1:nsteps_1D
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
for i=1:nsteps_1D
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

%if nsteps_1D > 1 %nsteps_1D gives a matlab error
%%Plot evolution of temperature profile
%figure;
%mytemp_i=double(mytemp_i);
%pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),mytemp_i(imin_1D:imax_1D,:))
%shading interp;
%colorbar;
%xlabel(timelab{1})
%ylabel(radlab)
%title('Evolution of normalized ion temperature')
%set(gca,'Units','centimeters','Position',[ax ay aw ah])
%end
%
%Plot evolution of R/LT
%
if adhoc == 1 & nsel_profile >1
RoLT=-dTdrhooa./mytemp_i*aspect_ratio;
else
RoLT=-dTds./mytemp_i*aspect_ratio;
end
RoLT=double(RoLT);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLT(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L_{T,i}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

f_av_i=double(f_av_i);
%if nsteps_1D > 1 %nsteps_1D gives a matlab error
%figure;
%pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),f_av_i(imin_1D:imax_1D,:))
%shading interp;
%colorbar;
%xlabel(timelab{1})
%ylabel(radlab)
%title('Evolution of normalized ion density')
%set(gca,'Units','centimeters','Position',[ax ay aw ah])
%end
%
%Plot evolution of R/Ln
%
if adhoc == 1 & nsel_profile >1
RoLn=-dndrhooa./f_av_i*aspect_ratio;
else
RoLn=-dnds./f_av_i*aspect_ratio;
end
RoLn=double(RoLn);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLn(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L{n,i}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAI IONS%%%%%%%%%%%%%%%%%%%%%
etai=RoLT./RoLn;
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),etai(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \eta_i')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM IONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing ion parallel momentum')
for i=1:nsteps_1D
int_v_par_av=csaps(double(psi_prof1D),double(v_par_av(:,i)),1-tol);
%v_par_av on a higher resolution \tilde(psi) grid
v_par_av_hr(:,i)=ppval(psi_prof1D_hr,int_v_par_av);
end

v_par_av_hr=double(v_par_av_hr);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr,v_par_av_hr)
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of ion parallel momentum')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing ion heatflux')
for i=1:nsteps_1D
int_efluxw_rad=csaps(double(psi_prof1D),double(efluxw_rad(:,i)),1-tol);
efluxw_rad_i(:,i)=ppval(psi_prof1D_hr,int_efluxw_rad);
end
%
%Plot evolution of Q
%
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
efluxw_rad_i=double(efluxw_rad_i);
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),efluxw_rad_i(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of ion normalized heat flux')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%initialize chi_perp
chi_perp=0*efluxw_rad_i;
if adhoc == 1 & nsel_profile >1 & norm_Dimits==1
chi_perp(imin_1D:imax_1D,:)=-efluxw_rad_i(imin_1D:imax_1D,:)./(f_av_i(imin_1D:imax_1D,:).*dTdrhooa(imin_1D:imax_1D,:))*(lx/2)^2/kappan0;
else
chi_perp(imin_1D:imax_1D,:)=-efluxw_rad_i(imin_1D:imax_1D,:)./(f_av_i(imin_1D:imax_1D,:).*dTds(imin_1D:imax_1D,:))*(lx/2)^2;
end
chi_perp=double(chi_perp);
if nsteps_1D > 1 %nsteps_1D=1 gives a matlab error
figure
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),chi_perp(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \chi_i')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
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

if nlelec == 0
figure;
hold on;
plot(RoLT_av_loc,chi_perp_av_loc,'kx-');
plot(RoLT_graph,chi_perp_Dimits,'b');
axis([RoLT_graph(1) RoLT_graph(end) 0 1.2*max(abs(chi_perp_av_loc))])
xlabel('R_0/L_{T,i}')
ylabel('\chi_i/\chi_{GB}')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE ELECTRONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlelec > 0
mytemp_el=melec/mion/3*(v_perp2_av_el./f_av_el+v_par2_av_el./f_av_el-(v_par_av_el./f_av_el).^2);
%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, electrons')
for i=1:nsteps_1D
int_mytemp_el=csaps(double(psi_prof1D),double(mytemp_el(:,i)),1-tol);
%coefficient of temperature gradient
dint_mytemp_el=fnder(int_mytemp_el);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i_el=ppval(psi_prof1D_hr,dint_mytemp_el);
%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i_el(:,i)=ppval(psi_prof1D_hr,int_mytemp_el);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds_el(:,i)=dmytemp_i_el.*s_prof1D_hr*2;
end

%
% Get dT/d(rho/a) on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=1:nsteps_1D
dTdrhooa_el(:,i)=dTds_el(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed temperature
figure;
plot(s_prof,te,'b',s_prof1D,mytemp_el(:,1),'r',s_prof1D_hr,mytemp_i_el(:,1),'k',s_prof1D,mytemp_el(:,nsteps_1D),'g',s_prof1D_hr,mytemp_i_el(:,nsteps_1D),'c')
xlabel('s')
ylabel('T_e/T_e(s_0)')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Initial electron temperature profile')

figure;
if adhoc == 1 & nsel_profile >1
%plot initial analytic and reconstructed logarithmic temperature gradient
plot(s_prof,gradte,'b',s_prof1D_hr,dTdrhooa_el(:,1)./mytemp_i_el(:,1),'r',s_prof1D_hr,dTdrhooa_el(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'k')
ylabel('a/L_{T,e,\rho/a}')
else
plot(s_prof,gradte,'b',s_prof1D_hr,dTds_el(:,1)./mytemp_i_el(:,1),'r',s_prof1D_hr,dTds_el(:,nsteps_1D)./mytemp_i_el(:,nsteps_1D),'k')
ylabel('a/L_{T,e,s}')
end
xlabel('s')
legend('analytic','particles, initial','particle, final')
title('Electron temperature gradient profile')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdatab=hdf5read(filename,'/equil/arrays/lambdatab');
alphatab=cos(lambdatab);
alphatab=max(alphatab,0);

sgrid_eq=hdf5read(filename,'/equil/arrays/sgrid_eq');
thgrid_eq=hdf5read(filename,'/equil/arrays/thgrid_eq');
nsgrid_eq=length(sgrid_eq);
nthgrid_eq=length(thgrid_eq);
dsgrid_eq=sgrid_eq(2)-sgrid_eq(1);
dthgrid_eq=thgrid_eq(2)-thgrid_eq(1);

%reconstruct equilibrium passing density
for i=1:length(s_prof)
is=find(sgrid_eq <= s_prof(i),1,'last');
if isempty(is)==1
is=1;
end
is=min(is,length(sgrid_eq)-1);
ws1=(s_prof(i)-sgrid_eq(is))/dsgrid_eq;
ws0=1-ws1;
alphab_fsa_int=ws0*alphab_fsa(is)+ws1*alphab_fsa(is+1);
ne0p(i)=ne(i)*(1-alphab_fsa_int);
ne0t(i)=ne(i)*alphab_fsa_int;
gradne0p(i)=gradne(i)*(1.0-alphab_fsa_int);
gradne0t(i)=gradne(i)*alphab_fsa_int;
end

%Get ne0p and gradne0p on s_prof1D grid
for i=1:length(s_prof1D)
is=find(s_prof <= s_prof1D(i),1,'last');
if isempty(is)==1
is=1;
end
is=min(is,length(s_prof)-1);
ws1=(s_prof1D(i)-s_prof(is))/(s_prof(is+1)-s_prof(is));
ws0=1-ws1;
ne0p_1D(i)=ws0*ne0p(is)+ws1*ne0p(is+1);
gradne0p_1D(i)=ws0*gradne0p(is)+ws1*gradne0p(is+1);
end
if size(ne0p_1D,2)>size(ne0p_1D,1)
ne0p_1D=transpose(ne0p_1D);
gradne0p_1D=transpose(gradne0p_1D);
end

if nlelec==2
for i=1:nsteps_1D
f_av_el_tot(:,i)=f_av_el(:,i)+ne0p_1D;
end
else
f_av_el_tot=f_av_el;
end

disp('Smoothing electron density and gradient')
for i=1:nsteps_1D
int_f_av_el=csaps(double(psi_prof1D),double(f_av_el_tot(:,i)),1-tol);
f_av_i_el_tot(:,i)=ppval(psi_prof1D_hr,int_f_av_el);
%coefficient of density gradient
dint_f_av_el=fnder(int_f_av_el);
%dn/d\tilde(psi) on a higher resolution \tilde(psi) grid
df_av_i_el=ppval(psi_prof1D_hr,dint_f_av_el);
%dn/ds, normalised on a higher resolution \tilde(psi) grid
dnds_el(:,i)=df_av_i_el.*s_prof1D_hr*2;         
end
%
% Get dn/d(rho/a) on a higher resolution \tilde(psi) grid
%
if adhoc ==1 & nsel_profile > 1
for i=1:nsteps_1D
dndrhooa_el(:,i)=dnds_el(:,i)./transpose(drhods_int);
end
end

%plot analytic and reconstructed density
figure;
plot(s_prof,ne,'b',s_prof1D,f_av_el_tot(:,1),'r',s_prof1D_hr,f_av_i_el_tot(:,1),'k',s_prof1D,f_av_el_tot(:,nsteps_1D),'g',s_prof1D_hr,f_av_i_el_tot(:,nsteps_1D),'c')
xlabel('s')
ylabel('n_e/<n_e>')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Electron density profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if adhoc == 1 & nsel_profile >1
plot(s_prof,gradne,'b',s_prof1D_hr,dndrhooa_el(:,1)./f_av_i_el_tot(:,1),'r',s_prof1D_hr,dndrhooa_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'k')
ylabel('a/L_{n,e,\rho/a}')
else
plot(s_prof,gradne,'b',s_prof1D_hr,dnds_el(:,1)./f_av_i_el_tot(:,1),'r',s_prof1D_hr, dnds_el(:,nsteps_1D)./f_av_i_el_tot(:,nsteps_1D),'k')
ylabel('a/L_{n,e,s}')
end
xlabel('s')
legend('analytic','particles, initial','particles, final')
title('Electron density gradient profile')

%if nsteps_1D > 1 %nsteps_1D gives a matlab error
%%Plot evolution of temperature profile
%figure;
%mytemp_i_el=double(mytemp_i_el);
%pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),mytemp_i_el(imin_1D:imax_1D,:))
%shading interp;
%colorbar;
%xlabel(timelab{1})
%ylabel(radlab)
%title('Evolution of normalized electron temperature')
%set(gca,'Units','centimeters','Position',[ax ay aw ah])
%end
%
%Plot evolution of R/LT
%
if adhoc == 1 & nsel_profile >1
RoLT_el=-dTdrhooa_el./mytemp_i_el*aspect_ratio;
else
RoLT_el=-dTds_el./mytemp_i_el*aspect_ratio;
end
RoLT_el=double(RoLT_el);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLT_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L_{T,e}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

f_av_i_el_tot=double(f_av_i_el_tot);
%if nsteps_1D > 1 %nsteps_1D gives a matlab error
%figure;
%pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),f_av_i_el(imin_1D:imax_1D,:))
%shading interp;
%colorbar;
%xlabel(timelab{1})
%ylabel(radlab)
%title('Evolution of normalized electron density')
%set(gca,'Units','centimeters','Position',[ax ay aw ah])
%end
%
%Plot evolution of R/Ln
%
if adhoc == 1 & nsel_profile >1
RoLn_el=-dndrhooa_el./f_av_i_el_tot*aspect_ratio;
else
RoLn_el=-dnds_el./f_av_i_el_tot*aspect_ratio;
end
RoLn_el=double(RoLn_el);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLn_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L{n,e}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM ELECTRONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing electron parallel momentum')
for i=1:nsteps_1D
int_v_par_av_el=csaps(double(psi_prof1D),double(v_par_av_el(:,i)),1-tol);
%v_par_av on a higher resolution \tilde(psi) grid
v_par_av_hr_el(:,i)=ppval(psi_prof1D_hr,int_v_par_av_el);
end

v_par_av_hr_el=double(v_par_av_hr_el);
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr,v_par_av_hr_el)
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron parallel momentum')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAE ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etae=RoLT_el./RoLn_el;
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),etae(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \eta_e')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing electron heatflux')
for i=1:nsteps_1D
int_efluxw_rad_el=csaps(double(psi_prof1D),double(efluxw_rad_el(:,i)),1-tol);
efluxw_rad_i_el(:,i)=ppval(psi_prof1D_hr,int_efluxw_rad_el);
end
%
%Plot evolution of Q
%
if nsteps_1D > 1 %nsteps_1D gives a matlab error
figure;
efluxw_rad_i_el=double(efluxw_rad_i_el);
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),efluxw_rad_i_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron normalized heat flux')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%initialize chi_perp
chi_perp_el=0*efluxw_rad_i_el;
if adhoc == 1 & nsel_profile >1 & norm_Dimits==1
chi_perp_el(imin_1D:imax_1D,:)=-efluxw_rad_i_el(imin_1D:imax_1D,:)./(f_av_i_el_tot(imin_1D:imax_1D,:).*dTdrhooa_el(imin_1D:imax_1D,:))*(lx/2)^2/kappan0;
else
chi_perp_el(imin_1D:imax_1D,:)=-efluxw_rad_i_el(imin_1D:imax_1D,:)./(f_av_i_el_tot(imin_1D:imax_1D,:).*dTds_el(imin_1D:imax_1D,:))*(lx/2)^2;
end
chi_perp_el=double(chi_perp_el);
if nsteps_1D > 1 %nsteps_1D=1 gives a matlab error
figure
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),chi_perp_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \chi_e')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%Plot chi vs R/LT
ilow=find(rad_prof1D_hr < radlow,1,'last');
iup =find(rad_prof1D_hr < radup,1,'last');
text=sprintf('%s', 'Profiles are averaged between ',radlab,'=');
disp(text);
rad_prof1D_hr(ilow)
rad_prof1D_hr(iup)

RoLT_av_el=squeeze(mean(RoLT_el(ilow:iup,:)));
chi_perp_av_el=squeeze(mean(chi_perp_el(ilow:iup,:)));
RoLT_av_loc_el(1)=RoLT_av_el(1);
chi_perp_av_loc_el(1)=chi_perp_av_el(1);

npoints=floor(nsteps_1D/tlength)-1;
for i=1:npoints
RoLT_av_loc_el(i+1)=RoLT_av_el(1+i*tlength);
chi_perp_av_loc_el(i+1)=chi_perp_av_el(1+i*tlength);
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
plot(RoLT_av_loc,chi_perp_av_loc,'bx-');
plot(RoLT_graph_el,chi_perp_Dimits_el,'k');
plot(RoLT_graph,chi_perp_Dimits,'b');
legend('electrons','ions')
xlabel('R_0/L_{T}')
ylabel('\chi/\chi_{GB}')
chimaxi=max(chi_perp_av_loc);
chimaxe=max(chi_perp_av_loc_el);
chimax=max(chimaxi, chimaxe);
RoLTmini=min(RoLT_av_loc);
RoLTmine=min(RoLT_av_loc_el);
RoLTmin=min(RoLTmini, RoLTmine);
RoLTmaxi=max(RoLT_av_loc);
RoLTmaxe=max(RoLT_av_loc_el);
RoLTmax=max(RoLTmaxi, RoLTmaxe);
axis([0.8*RoLTmin 1.2*RoLTmax 0 1.2*chimax])
end

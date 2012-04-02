function nemorb_transport_hdf5(sim, ind)
%################################
%nemorb_transport_transport(sim, ind)
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

if(~(size(sim(k).f_av,2)==length(sim(k).time_1D)))
sim(k).f_av=sim(k).f_av';
end

if(~(size(sim(k).v_par_av,2)==length(sim(k).time_1D)))
sim(k).v_par_av=sim(k).v_par_av';
end

if(~(size(sim(k).v_perp2_av,2)==length(sim(k).time_1D)))
sim(k).v_perp2_av=sim(k).v_perp2_av';
end

if(~(size(sim(k).v_par2_av,2)==length(sim(k).time_1D)))
sim(k).v_par2_av=sim(k).v_par2_av';
end

if sim(k).nsel_trap_only > 0
sim(k).f_av_el=sim(k).f_av_el + sim(k).v_p_el+sim(k).f_av_tr_el;
sim(k).v_par_av_el=sim(k).v_par_asim(k).v_p_el+sim(k).v_par_av_tr_el;
sim(k).v_perp2_av_el=sim(k).v_perp2_asim(k).v_p_el+sim(k).v_perp2_av_tr_el;
sim(k).v_par2_av_el=sim(k).v_par2_asim(k).v_p_el+sim(k).v_par2_av_tr_el;

if(~(size(sim(k).f_av_el,2)==length(sim(k).time_1D)))
sim(k).f_av_el=sim(k).f_av_el';
end

if(~(size(sim(k).v_par_av_el,2)==length(sim(k).time_1D)))
sim(k).v_par_av_el=sim(k).v_par_av_el';
end

if(~(size(sim(k).v_perp2_av_el,2)==length(sim(k).time_1D)))
sim(k).v_perp2_av_el=sim(k).v_perp2_av_el';
end

if(~(size(sim(k).v_par2_av_el,2)==length(sim(k).time_1D)))
sim(k).v_par2_av_el=sim(k).v_par2_av_el';
end
end

% $$$ if(~(size(sim(k).efluxw_rad,2)==length(sim(k).time_1D)))
% $$$ sim(k).efluxw_rad=sim(k).efluxw_rad';
% $$$ end

if sim(k).nsel_trap_only > 0
sim(k).efluxw_rad_el=sim(k).efluxw_rad_tr_el+sim(k).efluxw_rad_p_el;
if(~(size(sim(k).efluxw_rad_el,2)==length(sim(k).time_1D)))
sim(k).efluxw_rad_el=sim(k).efluxw_rad_el';
end
end

mytemp=1./3*(sim(k).v_perp2_av./sim(k).f_av+sim(k).v_par2_av./sim(k).f_av-(sim(k).v_par_av./sim(k).f_av).^2);
psifmin=sim(k).sfmin*sim(k).sfmin;
psifmax=sim(k).sfmax*sim(k).sfmax;
sim(k).psi_prof1D_hr=[psifmin:(psifmax-psifmin)/(sim(k).nhr-1):psifmax];
if(abs(sim(k).psi_prof1D_hr(end)-psifmax)>1e-6)
  sim(k).psi_prof1D_hr=[sim(k).psi_prof1D_hr psifmax];
end
sim(k).s_prof1D=sqrt(sim(k).psi_prof1D);
sim(k).s_prof1D_hr=sqrt(sim(k).psi_prof1D_hr);
sim(k).s_prof1D_hr=double(sim(k).s_prof1D_hr);

sim(k).nprof_1D=size(sim(k).v_perp2_av,1);
sim(k).nsteps_1D=size(sim(k).v_perp2_av,2);

%
% Get d(rho/a)/ds on a higher resolution \tilde(psi) grid
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile > 1
for i=1:sim(k).nhr
%get rho/a
psi_index=find(sim(k).psi_prof<sim(k).psi_prof1D_hr(i),1,'last');
if isempty(psi_index) 
psi_index=1;
end
w1=(sim(k).psi_prof1D_hr(i)-sim(k).psi_prof(psi_index))/(sim(k).psi_prof(psi_index+1)-sim(k).psi_prof(psi_index));
w0=1-w1;
%rho_int=rho/a (sim(k).psi_prof1D_hr(i))
rho_int=w0*sim(k).rho_prof(psi_index)+w1*sim(k).rho_prof(psi_index+1);
vol_prof1D_hr(i)=rho_int;
rho_index=find(sim(k).rho_prof< rho_int,1,'last');
if isempty(rho_index)
rho_index=1;
end
w1=(rho_int-sim(k).rho_prof(rho_index))/(sim(k).rho_prof(rho_index+1)-sim(k).rho_prof(rho_index));
w0=1-w1;
sim(k).drhods_int(i)=w0*sim(k).drhods(rho_index)+w1*sim(k).drhods(rho_index+1);
end
end

%select radial coordinates for plot
switch sim(k).irad_coord
case 1
rad_prof1D_hr=sim(k).s_prof1D_hr;
radlab='s';
case 2
rad_prof1D_hr=vol_prof1D_hr;
radlab='sqrt(V(s)/V)';
end

imin_1D=find(rad_prof1D_hr < sim(k).radmin_1D,1,'last');
imax_1D=find(rad_prof1D_hr < sim(k).radmax_1D,1,'last');

%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, ions')
for i=1:sim(k).nsteps_1D
%size(psi_prof1D) = 6, size(mytemp) = 66 3000)
int_mytemp=csaps(double(sim(k).psi_prof1D),double(mytemp(:,i)'),1-sim(k).tol);
%coefficient of temperature gradient
dint_mytemp=fnder(int_mytemp);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i=ppval(sim(k).psi_prof1D_hr,dint_mytemp);
%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i(:,i)=ppval(sim(k).psi_prof1D_hr,int_mytemp);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds(:,i)=dmytemp_i.*sim(k).s_prof1D_hr*2;
end

disp('done smoothing of T and grad_T')

%
% Get dT(rho/a)/d(rho/a) on a higher resolution \tilde(psi) grid
%

if sim(k).nsel_equil == 2 & sim(k).nsel_profile > 1
for i=1:sim(k).nsteps_1D
dTdrhooa(:,i)=dTds(:,i)./transpose(sim(k).drhods_int);
end
end

%plot analytic and reconstructed temperature
figure;
plot(sim(k).s_prof,sim(k).t_pic_MI,'b',sim(k).s_prof1D,mytemp(:,1),'r',sim(k).s_prof1D_hr,mytemp_i(:,1),'k',sim(k).s_prof1D,mytemp(:,sim(k).nsteps_1D),'g',sim(k).s_prof1D_hr,mytemp_i(:,sim(k).nsteps_1D),'c');
xlabel('s')
ylabel('T_i/T_e(s_0)')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Initial ion temperature profile')

figure;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
%plot initial analytic and reconstructed logarithmic temperature gradient
plot(sim(k).s_prof,sim(k).gradt_pic_MI,'b',sim(k).s_prof1D_hr,dTdrhooa(:,1)./mytemp_i(:,1),'r',sim(k).s_prof1D_hr,dTdrhooa(:,sim(k).nsteps_1D)./mytemp_i(:,sim(k).nsteps_1D),'k')
ylabel('a/L_{T,i,\rho/a}')
else
plot(sim(k).s_prof,sim(k).gradt_pic_MI,'b',sim(k).s_prof1D_hr,dTds(:,2)./mytemp_i(:,2),'r',sim(k).s_prof1D_hr,dTds(:,sim(k).nsteps_1D)./mytemp_i(:,sim(k).nsteps_1D),'k')
ylabel('a/L_{T,i,s}')
end
xlabel('s')
legend('analytic','particles, initial','particle, final')
title('Ion temperature gradient profile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Smoothing density and gradient')
  for i=1:sim(k).nsteps_1D
intf_av=csaps(double(sim(k).psi_prof1D),double(sim(k).f_av(:,i)),1-sim(k).tol);
sim(k).f_av_i(:,i)=ppval(sim(k).psi_prof1D_hr,intf_av);
%coefficient of density gradient
dintf_av=fnder(intf_av);
%dn/d\tilde(psi) on a higher resolution \tilde(psi) grid
dsim(k).f_av_i=ppval(sim(k).psi_prof1D_hr,dintf_av);
%dn/ds, normalised on a higher resolution \tilde(psi) grid
dnds(:,i)=dsim(k).f_av_i.*sim(k).s_prof1D_hr*2;         
end
%
% Get dn/d(rho/a) on a higher resolution \tilde(psi) grid
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile > 1
for i=1:sim(k).nsteps_1D
dndrhooa(:,i)=dnds(:,i)./transpose(sim(k).drhods_int);
end
end

%plot analytic and reconstructed density
figure;
plot(sim(k).s_prof,sim(k).n_pic_MI,'b',sim(k).s_prof1D,sim(k).f_av(:,1),'r',sim(k).s_prof1D_hr,sim(k).f_av_i(:,1),'k',sim(k).s_prof1D,sim(k).f_av(:,sim(k).nsteps_1D),'g',sim(k).s_prof1D_hr,sim(k).f_av_i(:,sim(k).nsteps_1D),'c')
xlabel('s')
ylabel('n_i/<n_i>')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Ion density profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
plot(sim(k).s_prof,sim(k).gradn_pic_MI,'b',sim(k).s_prof1D_hr,dndrhooa(:,1)./sim(k).f_av_i(:,1),'r',sim(k).s_prof1D_hr,dndrhooa(:,sim(k).nsteps_1D)./sim(k).f_av_i(:,sim(k).nsteps_1D),'k');
ylabel('a/L_{n,i,\rho/a}')
else
plot(sim(k).s_prof,sim(k).gradn_pic_MI,'b',sim(k).s_prof1D_hr,dnds(:,1)./sim(k).f_av_i(:,1),'r',sim(k).s_prof1D_hr, dnds(:,sim(k).nsteps_1D)./sim(k).f_av_i(:,sim(k).nsteps_1D),'k');
ylabel('a/L_{n,i,s}')
end
xlabel('s')
legend('analytic','particles, initial','particles, final')
title('Ion density gradient profile')

%
%Plot evolution of R/LT
%
sim(k).aspect_ratio = sim(k).r0_mid/sim(k).a_mid;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
RoLT=-dTdrhooa./mytemp_i*sim(k).aspect_ratio;
else
RoLT=-dTds./mytemp_i*sim(k).aspect_ratio;
end

%
%Plot evolution of R/Ln
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
RoLn=-dndrhooa./sim(k).f_av_i*sim(k).aspect_ratio;
else
RoLn=-dnds./sim(k).f_av_i*sim(k).aspect_ratio;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAI IONS%%%%%%%%%%%%%%%%%%%%%
etai=RoLT./RoLn;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM IONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing ion parallel momentum')
for i=1:sim(k).nsteps_1D
intv_par_av=csaps(double(sim(k).psi_prof1D),double(sim(k).v_par_av(:,i)),1-sim(k).tol);
%sim(k).v_par_av on a higher resolution \tilde(psi) grid
sim(k).v_par_av_hr(:,i)=ppval(sim(k).psi_prof1D_hr,intv_par_av);
end

%sim(k).v_par_av_hr=sim(k).v_par_av_hr./sim(k).f_av_i;
%Now follows some guesswork from chi ions

chi_phi(imin_1D:imax_1D,:)=-sim(k).v_par_av_hr(imin_1D:imax_1D,:)./(sim(k).f_av_i(imin_1D:imax_1D,:).*dnds(imin_1D:imax_1D,:))*(sim(k).lx/2)^2;

figure;
plot(chi_phi(:,1))
hold on;
plot(chi_phi(:,end))
hold off;

%Plot chi vs R/Ln
ilow=find(rad_prof1D_hr < sim(k).radlow,1,'last');
iup =find(rad_prof1D_hr < sim(k).radup,1,'last');
text=sprintf('%s', 'Profiles are averaged between ',radlab,'=');
disp(text);
rad_prof1D_hr(ilow);
rad_prof1D_hr(iup);

RoLn_av=squeeze(mean(RoLn(ilow:iup,:)));
chi_phi_av=squeeze(mean(chi_phi(ilow:iup,:)));
RoLn_av_loc(1)=RoLn_av(1);
chi_phi_av_loc(1)=chi_phi_av(1);

npoints=floor(sim(k).nsteps_1D/sim(k).tlength)-1;
for i=1:npoints
RoLn_av_loc(i+1)=RoLn_av(1+i*sim(k).tlength);
chi_phi_av_loc(i+1)=chi_phi_av(1+i*sim(k).tlength);
end

RoLn_graph=[0.8*min(RoLn_av_loc):0.01:max(4.1,1.2*max(RoLn_av_loc))];
if sim(k).norm_Dimits==1
chi_phi_Dimits=15.4*(1-6./RoLn_graph);
else
chi_phi_Dimits=15.4*sim(k).kappan0*(1-6./RoLn_graph);
end

if sim(k).nsel_trap_only == 0
figure;
if(~exist('uptotime')) uptotime=1.0;end
hold on;
plot(RoLn_av_loc(1:round(end*uptotime)),chi_phi_av_loc(1:round(end*uptotime)),'kx-');
plot(RoLn_graph,chi_phi_Dimits,'b');
hold off
axis([RoLn_graph(1) RoLn_graph(end) 0 1.2*max(abs(chi_phi_av_loc))])
xlabel('R_0/L_{n,i}')
ylabel('\chi_{phi}/\chi_{GB}')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI IONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing ion heatflux')
for i=1:sim(k).nsteps_1D
intefluxw_rad=csaps(double(sim(k).psi_prof1D),double(sim(k).efluxw_rad(:,i)),1-sim(k).tol);
sim(k).efluxw_rad_i(:,i)=ppval(sim(k).psi_prof1D_hr,intefluxw_rad);
end
%initialize chi_perp

%chi_perp=0*sim(k).efluxw_rad_i;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1 & norm_Dimits==1
chi_perp(imin_1D:imax_1D,:)=-sim(k).efluxw_rad_i(imin_1D:imax_1D,:)./(sim(k).f_av_i(imin_1D:imax_1D,:).*dTdrhooa(imin_1D:imax_1D,:))*(sim.lx/2)^2/kappan0;
else
chi_perp(imin_1D:imax_1D,:)=-sim(k).efluxw_rad_i(imin_1D:imax_1D,:)./(sim(k).f_av_i(imin_1D:imax_1D,:).*dTds(imin_1D:imax_1D,:))*(sim.lx/2)^2;
end
%Plot chi vs R/LT
ilow=find(rad_prof1D_hr < sim(k).radlow,1,'last');
iup =find(rad_prof1D_hr < sim(k).radup,1,'last');
text=sprintf('%s', 'Profiles are averaged between ',radlab,'=');
disp(text);
rad_prof1D_hr(ilow)
rad_prof1D_hr(iup)

RoLT_av=squeeze(mean(RoLT(ilow:iup,:)));
chi_perp_av=squeeze(mean(chi_perp(ilow:iup,:)));
RoLT_av_loc(1)=RoLT_av(1);
chi_perp_av_loc(1)=chi_perp_av(1);

npoints=floor(sim(k).nsteps_1D/sim(k).tlength)-1;
for i=1:npoints
RoLT_av_loc(i+1)=RoLT_av(1+i*sim(k).tlength);
chi_perp_av_loc(i+1)=chi_perp_av(1+i*sim(k).tlength);
end

RoLT_graph=[0.8*min(RoLT_av_loc):0.01:max(4.1,1.2*max(RoLT_av_loc))];
if sim(k).norm_Dimits==1
chi_perp_Dimits=15.4*(1-6./RoLT_graph);
else
chi_perp_Dimits=15.4*sim(k).kappan0*(1-6./RoLT_graph);
end

if sim(k).nsel_trap_only == 0
figure;
if(~exist('uptotime')) uptotime=1.0;end
hold on;
plot(RoLT_av_loc(1:round(end*uptotime)),chi_perp_av_loc(1:round(end*uptotime)),'kx-');
plot(RoLT_graph,chi_perp_Dimits,'b');
hold off
axis([RoLT_graph(1) RoLT_graph(end) 0 1.2*max(abs(chi_perp_av_loc))])
xlabel('R_0/L_{T,i}')
ylabel('\chi_i/\chi_{GB}')
end

% Electrons below this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE ELECTRONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim(k).nsel_trap_only > 0
mytemp_el=melec/mion/3*(sim(k).v_perp2_av_el./sim(k).f_av_el+sim(k).v_par2_av_el./sim(k).f_av_el-(sim(k).v_par_av_el./sim(k).f_av_el).^2);
%
% Smoothing (cubic splines) of temperature and temperature gradient 
%
disp('Smoothing temperature and gradient, electrons')
for i=1:sim(k).nsteps_1D
int_mytemp_el=csaps(double(sim(k).psi_prof1D),double(mytemp_el(:,i)),1-sim(k).tol);
%coefficient of temperature gradient
dint_mytemp_el=fnder(int_mytemp_el);
%dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
dmytemp_i_el=ppval(sim(k).psi_prof1D_hr,dint_mytemp_el);
%smoothed temperature, normalised on a higher resolution \tilde(psi) grid
mytemp_i_el(:,i)=ppval(sim(k).psi_prof1D_hr,int_mytemp_el);
%dT/ds, normalised on a higher resolution \tilde(psi) grid
dTds_el(:,i)=dmytemp_i_el.*sim(k).s_prof1D_hr*2;
end

%
% Get dT/d(rho/a) on a higher resolution \tilde(psi) grid
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile > 1
for i=1:sim(k).nsteps_1D
dTdrhooa_el(:,i)=dTds_el(:,i)./transpose(sim(k).drhods_int);
end
end

%plot analytic and reconstructed temperature
figure;
plot(sim(k).s_prof,te,'b',sim(k).s_prof1D,mytemp_el(:,1),'r',sim(k).s_prof1D_hr,mytemp_i_el(:,1),'k',sim(k).s_prof1D,mytemp_el(:,sim(k).nsteps_1D),'g',sim(k).s_prof1D_hr,mytemp_i_el(:,sim(k).nsteps_1D),'c')
xlabel('s')
ylabel('T_e/T_e(s_0)')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Initial electron temperature profile')

figure;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
%plot initial analytic and reconstructed logarithmic temperature gradient
plot(sim(k).s_prof,gradte,'b',sim(k).s_prof1D_hr,dTdrhooa_el(:,1)./mytemp_i_el(:,1),'r',sim(k).s_prof1D_hr,dTdrhooa_el(:,sim(k).nsteps_1D)./mytemp_i_el(:,sim(k).nsteps_1D),'k')
ylabel('a/L_{T,e,\rho/a}')
else
plot(sim(k).s_prof,gradte,'b',sim(k).s_prof1D_hr,dTds_el(:,1)./mytemp_i_el(:,1),'r',sim(k).s_prof1D_hr,dTds_el(:,sim(k).nsteps_1D)./mytemp_i_el(:,sim(k).nsteps_1D),'k')
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
for i=1:length(sim(k).s_prof)
is=find(sgrid_eq <= sim(k).s_prof(i),1,'last');
if isempty(is)==1
is=1;
end
is=min(is,length(sgrid_eq)-1);
ws1=(sim(k).s_prof(i)-sgrid_eq(is))/dsgrid_eq;
ws0=1-ws1;
alphab_fsa_int=ws0*alphab_fsa(is)+ws1*alphab_fsa(is+1);
ne0p(i)=ne(i)*(1-alphab_fsa_int);
ne0t(i)=ne(i)*alphab_fsa_int;
gradne0p(i)=gradne(i)*(1.0-alphab_fsa_int);
gradne0t(i)=gradne(i)*alphab_fsa_int;
end

%Get ne0p and gradne0p on sim(k).s_prof1D grid
for i=1:length(sim(k).s_prof1D)
is=find(sim(k).s_prof <= sim(k).s_prof1D(i),1,'last');
if isempty(is)==1
is=1;
end
is=min(is,length(sim(k).s_prof)-1);
ws1=(sim(k).s_prof1D(i)-sim(k).s_prof(is))/(sim(k).s_prof(is+1)-sim(k).s_prof(is));
ws0=1-ws1;
ne0p_1D(i)=ws0*ne0p(is)+ws1*ne0p(is+1);
gradne0p_1D(i)=ws0*gradne0p(is)+ws1*gradne0p(is+1);
end
if size(ne0p_1D,2)>size(ne0p_1D,1)
ne0p_1D=transpose(ne0p_1D);
gradne0p_1D=transpose(gradne0p_1D);
end

if sim(k).nsel_trap_only==2
for i=1:sim(k).nsteps_1D
sim(k).f_av_el_tot(:,i)=sim(k).f_av_el(:,i)+ne0p_1D;
end
else
sim(k).f_av_el_tot=sim(k).f_av_el;
end

disp('Smoothing electron density and gradient')
for i=1:sim(k).nsteps_1D
intf_av_el=csaps(double(sim(k).psi_prof1D),double(sim(k).f_av_el_tot(:,i)),1-sim(k).tol);
sim(k).f_av_i_el_tot(:,i)=ppval(sim(k).psi_prof1D_hr,intf_av_el);
%coefficient of density gradient
dintf_av_el=fnder(intf_av_el);
%dn/d\tilde(psi) on a higher resolution \tilde(psi) grid
dsim(k).f_av_i_el=ppval(sim(k).psi_prof1D_hr,dintf_av_el);
%dn/ds, normalised on a higher resolution \tilde(psi) grid
dnds_el(:,i)=dsim(k).f_av_i_el.*sim(k).s_prof1D_hr*2;         
end
%
% Get dn/d(rho/a) on a higher resolution \tilde(psi) grid
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile > 1
for i=1:sim(k).nsteps_1D
dndrhooa_el(:,i)=dnds_el(:,i)./transpose(sim(k).drhods_int);
end
end

%plot analytic and reconstructed density
figure;
plot(sim(k).s_prof,ne,'b',sim(k).s_prof1D,sim(k).f_av_el_tot(:,1),'r',sim(k).s_prof1D_hr,sim(k).f_av_i_el_tot(:,1),'k',sim(k).s_prof1D,sim(k).f_av_el_tot(:,sim(k).nsteps_1D),'g',sim(k).s_prof1D_hr,sim(k).f_av_i_el_tot(:,sim(k).nsteps_1D),'c')
xlabel('s')
ylabel('n_e/<n_e>')
legend('analytic','particles, initial','smoothed, initial','particles, final','smoothed, final')
title('Electron density profile')

%plot analytic and reconstructed logarithmic density gradient
figure;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
plot(sim(k).s_prof,gradne,'b',sim(k).s_prof1D_hr,dndrhooa_el(:,1)./sim(k).f_av_i_el_tot(:,1),'r',sim(k).s_prof1D_hr,dndrhooa_el(:,sim(k).nsteps_1D)./sim(k).f_av_i_el_tot(:,sim(k).nsteps_1D),'k')
ylabel('a/L_{n,e,\rho/a}')
else
plot(sim(k).s_prof,gradne,'b',sim(k).s_prof1D_hr,dnds_el(:,1)./sim(k).f_av_i_el_tot(:,1),'r',sim(k).s_prof1D_hr, dnds_el(:,sim(k).nsteps_1D)./sim(k).f_av_i_el_tot(:,sim(k).nsteps_1D),'k')
ylabel('a/L_{n,e,s}')
end
xlabel('s')
legend('analytic','particles, initial','particles, final')
title('Electron density gradient profile')

%
%Plot evolution of R/LT
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
RoLT_el=-dTdrhooa_el./mytemp_i_el*sim(k).aspect_ratio;
else
RoLT_el=-dTds_el./mytemp_i_el*sim(k).aspect_ratio;
end
%
%Plot evolution of R/Ln
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
RoLn_el=-dndrhooa_el./sim(k).f_av_i_el_tot*sim(k).aspect_ratio;
else
RoLn_el=-dnds_el./sim(k).f_av_i_el_tot*sim(k).aspect_ratio;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL MOMENTUM ELECTRONS%%%%%%%%%%%%%%%%%%%%%
%
% Smoothing (cubic splines) of parallel momentum
%
disp('Smoothing electron parallel momentum')
for i=1:sim(k).nsteps_1D
intv_par_av_el=csaps(double(sim(k).psi_prof1D),double(sim(k).v_par_av_el(:,i)),1-sim(k).tol);
%sim(k).v_par_av on a higher resolution \tilde(psi) grid
sim(k).v_par_av_hr_el(:,i)=ppval(sim(k).psi_prof1D_hr,intv_par_av_el);
end

sim(k).v_par_av_hr_el=sim(k).v_par_av_hr_el./sim(k).f_av_i_el_tot;
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETAE ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etae=RoLT_el./RoLn_el;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHI ELECTRONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Heat Flux Smoothing
%
disp('Smoothing electron heatflux')
for i=1:sim(k).nsteps_1D
intefluxw_rad_el=csaps(double(sim(k).psi_prof1D),double(sim(k).efluxw_rad_el(:,i)),1-sim(k).tol);
sim(k).efluxw_rad_i_el(:,i)=ppval(sim(k).psi_prof1D_hr,intefluxw_rad_el);
end

%initialize chi_perp
chi_perp_el=0*sim(k).efluxw_rad_i_el;
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1 & norm_Dimits==1
chi_perp_el(imin_1D:imax_1D,:)=-sim(k).efluxw_rad_i_el(imin_1D:imax_1D,:)./(sim(k).f_av_i_el_tot(imin_1D:imax_1D,:).*dTdrhooa_el(imin_1D:imax_1D,:))*(sim.lx/2)^2/kappan0;
else
chi_perp_el(imin_1D:imax_1D,:)=-sim(k).efluxw_rad_i_el(imin_1D:imax_1D,:)./(sim(k).f_av_i_el_tot(imin_1D:imax_1D,:).*dTds_el(imin_1D:imax_1D,:))*(sim.lx/2)^2;
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

npoints=floor(sim(k).nsteps_1D/sim(k).tlength)-1;
for i=1:npoints
RoLT_av_loc_el(i+1)=RoLT_av_el(1+i*sim(k).tlength);
chi_perp_av_loc_el(i+1)=chi_perp_av_el(1+i*sim(k).tlength);
end

RoLT_graph_el=[0.8*min(RoLT_av_loc_el):0.01:max(4.1,1.2*max(RoLT_av_loc_el))];
if norm_Dimits==1
chi_perp_Dimits_el=15.4*(1-6./RoLT_graph_el);
else
chi_perp_Dimits_el=15.4*kappan0*(1-6./RoLT_graph_el);
end

% End of electrons

% Graphs - inc. electrons
figure;
hold on;
plot(RoLT_av_loc_el,chi_perp_av_loc_el,'rx-');
plot(RoLT_av_loc,chi_perp_av_loc,'bx-');
plot(RoLT_graph_el,chi_perp_Dimits_el,'k');
plot(RoLT_graph,chi_perp_Dimits,'k');
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


end % simulation loop
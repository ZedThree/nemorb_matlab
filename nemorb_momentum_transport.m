
function nemorb_momentum_transport(sim, ind)
%################################
%nemorb_momentum_transport(sim, ind)
%################################
%---------------
%
%Description
%
%Do the mometum transport and suchlike
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

rad_prof1D_hr = sim(k).s_prof1D_hr;
radlab = 's';

imin_1D=find(rad_prof1D_hr < sim(k).radmin_1D,1,'last');
imax_1D=find(rad_prof1D_hr < sim(k).radmax_1D,1,'last');

% $$$ %
% $$$ % Smoothing (cubic splines) of temperature and temperature gradient 
% $$$ %
% $$$ disp('Smoothing temperature and gradient, ions')
% $$$ for i=1:sim(k).nsteps_1D
% $$$ %size(psi_prof1D) = 6, size(mytemp) = 66 3000)
% $$$ int_mytemp=csaps(double(sim(k).psi_prof1D),double(mytemp(:,i)'),1-sim(k).tol);
% $$$ %coefficient of temperature gradient
% $$$ dint_mytemp=fnder(int_mytemp);
% $$$ %dT/d\tilde(psi) on a higher resolution \tilde(psi) grid
% $$$ dmytemp_i=ppval(sim(k).psi_prof1D_hr,dint_mytemp);
% $$$ %smoothed temperature, normalised on a higher resolution \tilde(psi) grid
% $$$ mytemp_i(:,i)=ppval(sim(k).psi_prof1D_hr,int_mytemp);
% $$$ %dT/ds, normalised on a higher resolution \tilde(psi) grid
% $$$ dTds(:,i)=dmytemp_i.*sim(k).s_prof1D_hr*2;
% $$$ end


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

% $$$ %
% $$$ %Plot evolution of R/LT
% $$$ %
% $$$ sim(k).aspect_ratio = sim(k).r0_mid/sim(k).a_mid;
% $$$ if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
% $$$ RoLT=-dTdrhooa./mytemp_i*sim(k).aspect_ratio;
% $$$ else
% $$$ RoLT=-dTds./mytemp_i*sim(k).aspect_ratio;
% $$$ end

%
%Plot evolution of R/Ln
%
if sim(k).nsel_equil == 2 & sim(k).nsel_profile >1
RoLn=-dndrhooa./sim(k).f_av_i*sim(k).aspect_ratio;
else
RoLn=-dnds./sim(k).f_av_i*sim(k).aspect_ratio;
end

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
xlabel('R_0/L_{T,i}')
ylabel('\chi_i/\chi_{GB}')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAPHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %Loop on simulations
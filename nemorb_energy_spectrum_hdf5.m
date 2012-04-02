function nemorb_energy_spectrum_hdf5(sim, ind)
%################################
%sortie=energy_spectrum_hdf5(sim,ind)
%################################
%---------------
%
%Description
%
%Plot energy of toroidal modes vs. time
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
%sortie = not used
%---------------

global timelab
if exist('ind') == 0
ind=1;
end

pwd_old=pwd;
for i=1:length(ind) %Loop on simulations
k=ind(i);
cd(sim(k).path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read poloidal and toroidal modes
n=hdf5read(sim(k).filename,'/data/var2d/generic/efspecmn1/coord1');
m=hdf5read(sim(k).filename,'/data/var2d/generic/efspecmn1/coord2');
nmodes=length(n);
mmodes=length(m);
%                              /data     /var2d/generic  /efspecmmn   /data       
time_ef=hdf5read(sim(k).filename,'/data/var2d/generic/efspecmn1/time');
nsteps_ef=length(time_ef);
efspecn=zeros(nmodes,nsteps_ef);
dt_ef=time_ef(2)-time_ef(1);

%Check that linear times for growth rate calculation are inside time interval
if sim(k).tmax > sim(k).time(end)
disp('tmax modified because it was too big')
sim(k).tmax=sim(k).time(end)
end

if sim(k).tmin > sim(k).time(end)
disp('tmin modified because it was too big')
sim(k).tmin=sim(k).time(1)
end
imin=1+floor(sim(k).tmin/dt_ef);
imax =1+floor(sim(k).tmax/dt_ef);

for j=1:length(n) %Loop on toroidal modes
nloc=n(j);
clear tmp;
%Get efspec(nloc,:,:)= (m,t) components of toroidal mode nloc
%warning in hdf5read_slice: dimensions in reverse order!!!
tmp=hdf5read_slice_new(sim(k).filename,'/data/var2d/generic/efspecmn1/data',[0 0 j-1],[nsteps_ef mmodes 1]);
%remove the 1st useless dimension
tmp=squeeze(tmp);
%tmp is a 2d array (m,t)
%Sum over poloidal modes to get the total toroidal energy
efspecn(j,:)=sum(tmp,1);
leg{j}=strcat('n=',num2str(nloc));
% $$$ Estart(j)=log10(efspecn(j,imin));
% $$$ Estop(j) =log10(efspecn(j,imax));
Estart(j)=(efspecn(j,imin));
Estop(j) =(efspecn(j,imax));
%Get linear growth rate
[coef, fit]=polyfit(double(time_ef(imin:imax)),transpose(log(efspecn(j,imin:imax))),1);
gamma(j)=coef(1)/2
end %Loop on toroidal modes

%Fastest growing mode
[gmax n0]=max(gamma);
n0max=n(n0);

%Plot energy of toroidal modes
figure;
if sim(k).nfilt1 == 0
plot(time_ef,log10(efspecn(1,:)),'k--');
plot(time_ef,log10(efspecn(2:end,:)))
else
semilogy(time_ef,efspecn(1:end,:))
hold on;
end
xlabel(timelab{1})
ylabel('E^{(n)} [m_ic_s^2]')
%title(sim(k).name)
%title(strcat('\gamma=',num2str(gamma(n0))));
title(sprintf('%s, \\gamma = %e',sim(k).name,gamma(n0)))
legend(leg,'Location','NorthEast')

t=[sim(k).tmin sim(k).tmax];
E=[Estart; Estop];
plot(t,E,'r:x');

%Plot growth rates of toroidal modes
% $$$ figure;
% $$$ plot(n,gamma,'bx-')
% $$$ xlabel('n')
% $$$ ylabel('\gamma [\Omega_i]')
end %loop on simulations
%sortie='OK';

cd(pwd_old);
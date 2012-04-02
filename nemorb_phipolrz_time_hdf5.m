function [potrz, potbarrz, fpot]=nemorb_phipolrz_time_hdf5(sim, t,lplot1, lplot2, ind)
%######################################
%[potrz, potbarrz]=phipolrz_time_hdf5(sim, t,lplot1, lplot2, ind)
%######################################
%---------------
%
%Description
%Plot phi and phi-phibar on poloidal plane
%---------------
%
%Dummy arguments
%
%sim   = array containing the simulations data, created with loader.m
%t     : time
%lplot1 : true = plot phi
%lplot2 : true = plot phibar
%ind   = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
%potrz = potential on the poloidal plane
%potbarrz = potential on the poloidal plane
%---------------

if exist('ind')==0
ind=1;
end

pwd_old=pwd;

potbarrz=0.0;
for i=1:length(ind) %loop on simulations
k=ind(i);
cd(sim(k).path)
time_potrz=hdf5read(sim(k).filename,'/data/var2d/generic/potrz/time');
nsteps_potrz=length(time_potrz);
%check if data for t=time exists;otherwise use nearest time
it=find(time_potrz <= t,1,'last');
treal=time_potrz(it);
if treal ~= t
disp('Warning: no data available for t=')
t
disp('Plot data for t=')
treal
end

%Read r and z coordinates
r=hdf5read(sim(k).filename,'/data/var2d/generic/potrz/coord1');
z=hdf5read(sim(k).filename,'/data/var2d/generic/potrz/coord2');
srz=hdf5read(sim(k).filename,'/data/var2d/generic/potrz/srz');
cont=[0:0.1:1];
nr=length(r);
nz=length(z);
nreq=size(srz,1);
nzeq=size(srz,2);
requ=[r(1):(r(end)-r(1))/(nreq-1):r(end)];
zequ=[z(1):(z(end)-z(1))/(nzeq-1):z(end)];
drequ=r(2)-r(1);
dzequ=z(2)-z(1);

%read data
potrz(1:nr,1:nz,i)=hdf5read_slice(sim(k).filename,'/data/var2d/generic/potrz/data',[it-1 0 0],[1 nr nz]);

%plot phi on (r,z) plane
if lplot1
figure;
hold on;
pcolor(double(r),double(z),double(transpose(potrz(1:nr,1:nz,i))));
shading interp;
%store color axis values
col=caxis;
%plot magnetic surfaces
contour(requ,zequ,srz',cont,'k--');
%restore old axis
caxis(col);colorbar;
xlabel('r [\rho_s]');
ylabel('z [\rho_s]');
axis equal image;
titlen=strcat('\phi at t=',num2str(treal),', ',sim(k).name);
title(titlen)
fpot(i)=gcf;
end

if sim(k).generic.nsel_phibar_solver ~= 99 %if n=0 is solved
%compute (phi-phibar) on (r,z) plane
time_phib_tmp=hdf5read(sim(k).filename,'/data/var1d/phibar/time');
%check if data exists for t=treal; otherwise take nearest time
t_phib=find(time)phib_tmp==treal);
if isempty(t_phib)==1
disp('cannot plot phi-\bar{phi} for t=')
treal
t_phib=find(time_phib_tmp<=treal,1,'last');
disp('compute phibar for t=')
time_phib_tmp(t_phib)
end

potbarrz=0*potrz;
%read data
sgrid=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/coord1');
phibar=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/data');
dsgrid=sgrid(2)-sgrid(1);

%get phibar(s) on the (r,z) grid by linear interpolation
potbarrz=potrz;
for j=1:size(potrz,1)
      ir=find(requ <= r(j),1,'last');
      ir=max(1,min(size(requ,2)-1,ir));
      wr1=(r(j)-requ(ir))/drequ;
      wr0=1.0-wr1;
   for m=1:size(potrz,2)
      iz=find(zequ <= z(m),1,'last');
      iz=max(1,min(size(zequ,2)-1,iz));
      wz1=(z(m)-zequ(iz))/dzequ;
      wz0=1.0-wz1;
      sloc=wr0*(wz0*srz(ir,iz)+wz1*srz(ir,iz+1))+wr1*(wz0*srz(ir+1,iz)+wz1*srz(ir+1,iz+1));
      sloc=max(0,min(1,sloc));
      is=find(sgrid <= sloc,1,'last');
      is=max(1,min(size(sgrid,1)-1,is));
      if isempty(is) == 1
      is=1;
      phibar_int=0.0;
      else
      ws1=(sloc-sgrid(is))/dsgrid;
      ws0=1.0-ws1;
      phibar_int=ws0*phibar(is,t_phib)+ws1*phibar(is+1,t_phib);
      end
      if is == 1 || is == 2
         potbarrz(j,m,i)=0;
      else
         potbarrz(j,m,i)=potbarrz(j,m,i)-phibar_int;
      end
    end
end

%plot (phi-phibar) on (r,z) plane
if lplot2
figure;
hold on;
pcolor(double(r),double(z),squeeze(double(transpose(potbarrz(:,:,i)))));
shading interp;
%store color axis values
col=caxis;
%plot magnetic surfaces
contour(requ,zequ,srz,cont,'k--');
%restore old axis
caxis(col);colorbar;
xlabel('r [\rho_s]');
ylabel('z [\rho_s]');
axis equal image;
titlen=strcat('\phi-bar{\phi} at t=',num2str(treal),', ',sim(k).name);
title(titlen)
end

end 

%output
tmp=double(squeeze(potrz(:,:,i)));
potrz(:,:,i)=tmp;

if sim(k).generic.nsel_phibar_solver ~=99
tmp=double(transpose(squeeze(potbarrz(:,:,i))));
potbarrz(:,:,i)=tmp;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%create a subplot with all the figure
%c=figure;
%nsim=length(ind);
%ncol=floor(sqrt(nsim));
%nrow=ceil(nsim/ncol);
%if nrow*ncol < nsim
%disp('error')
%else
%
%for i=1:length(ind)
%b(i)=subplot(nrow,ncol,i)
%baxis(i)=gca
%bpos(1:4,i)=get(baxis(i),'Position');
%end
%close(c);
%for i=1:length(ind)
%figure(fpot(i))
%faxis(i)=gca;
%set(faxis(i),'Position',bpos(1:4,i));
%end
%clear c
%c=figure;
%for i=1:length(ind)
%set(faxis(i),'parent',c)
%end
%end
cd(pwd_old);

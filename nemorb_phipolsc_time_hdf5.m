function [potsc, potbarsc, fpot]=nemorb_phipolsc_time_hdf5(sim, t,lplot1, lplot2, ind)
%######################################
%[potsc, potbarsc]=phipolsc_time_hdf5(sim, t,lplot1, lplot2, ind)
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
%potsc = potential on the poloidal plane
%potbarsc = potential on the poloidal plane
%---------------

if exist('ind')==0
  ind=1;
end

pwd_old=pwd;

potbarsc=0.0;
for i=1:length(ind) %loop on simulations
  k=ind(i);
  cd(sim(k).path)
  time_potsc=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/time');
  nsteps_potsc=length(time_potsc);
  %check if data for t=time exists;otherwise use nearest time
  it=find(time_potsc <= t,1,'last');
  treal=time_potsc(it);
  if treal ~= t
    disp('Warning: no data available for t=')
    t
    disp('Plot data for t=')
    treal
  end

  %Read r and z coordinates
  r=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/rsc');
  z=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/zsc');
  s=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/coord1');
  chi=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/coord2');
  cont=[0:0.1:1];
  ns=length(s);
  nchi=length(chi);
  s2=repmat(s,[1 nchi]);
% $$$ nreq=size(srz,1);
% $$$ nzeq=size(srz,2);
% $$$ requ=[r(1):(r(end)-r(1))/(nreq-1):r(end)];
% $$$ zequ=[z(1):(z(end)-z(1))/(nzeq-1):z(end)];
% $$$ drequ=r(2)-r(1);
% $$$ dzequ=z(2)-z(1);

  %read data
  potsc(1:ns,1:nchi,k)=hdf5read_slice_new(sim(k).filename,'/data/var2d/generic/potsc/data',[it-1 0 0],[1 nchi ns]);

  %plot phi on (r,z) plane
  if lplot1
    figure('Name',[sim(k).name ' - Phipol plot']);
    hold on;
    pcolor(r, z, potsc(:,:,k));
    shading interp;
    %store color axis values
    col=caxis;
    %plot magnetic surfaces
    contour(r,z,s2,cont,'k--');
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
    time_1D=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/time');
    %check if data exists for t=treal; otherwise take nearest time
    t_phib=find(time_1D==treal);
    if isempty(t_phib)==1
      disp('cannot plot phi-\bar{phi} for t=')
      treal
      t_phib=find(time_1D<=treal,1,'last');
      disp('compute phibar for t=')
      time_1D(t_phib)
    end

    potbarsc=0*potsc;
    %read data
    sgrid=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/coord1');
    phibar=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/data');
    %phibar_back=hdf5read(sim(k).filename,'/data/var1d/generic/phibar/phibar_back');
    %phibar=phibar-(repmat(phibar_back,1,size(phibar,2)));
    dsgrid=sgrid(2)-sgrid(1);
    phibar_hires(:,:,k)=interp1(sgrid,phibar(:,t_phib),s);
    potbarsc(:,:,k)=repmat(phibar_hires(:,:,k),[1 nchi]);

    %plot (phi-phibar) on (r,z) plane
    if lplot2
      figure('Name',[sim(k).name ' - Phi-phibar plot']);
      hold on;
      pcolor(r,z,double(potsc(:,:,k)-potbarsc(:,:,k)));
      shading interp;
      %store color axis values
      col=caxis;
      %plot magnetic surfaces
      contour(r,z,s2,cont,'k--');
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
% $$$ tmp=double((squeeze(potsc(:,:,i))));
% $$$ potsc(:,:,i)=tmp;
% $$$ 
% $$$ if sim(k).generic.nsel_phibar_solver ~=99
% $$$ tmp=double((squeeze(potbarsc(:,:,i))));
% $$$ potbarsc(:,:,i)=tmp;
% $$$ end

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

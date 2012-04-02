function sortie=plot_spec_hdf5(sim,species,ispec,iloc,n0,t_spec)
global qmin qmax qs0 deltam;

filename=sim.filename;
cd(sim.path);
%select which spectrum to read
switch ispec
case 1
varname='efspec1mn';
varlabel='E';
case 2
varname='phispecmn';
varlabel='\phi^2';
case 3
varname='rhospecmn';
varlabel='\rho^2';
end

%choose global/local spectrum
if iloc==1
varname=strcat(varname,'_loc');
typename='Local ';
maxis_min=-n0*qs0-deltam;
maxis_max=-n0*qs0+deltam;
else
typename='Global ';
qmax=sim.generic.q(end);
qmin=sim.generic.q(1);
deltam=10;
maxis_min=-n0*qmax-deltam
maxis_max=-n0*qmin+deltam
end

%define pathname
pathname=strcat('/data/var2d/',species,'/',varname,'/data');

%read (m,n) modes of spectrum
n=hdf5read(filename,strcat('/data/var2d/',species,'/',varname,'/coord1'));
m=hdf5read(filename,strcat('/data/var2d/',species,'/',varname,'/coord2'));
mmodes=length(m);
nmodes=length(n);

%read time
time_spec=hdf5read(filename,strcat('/data/var2d/',species,'/',varname,'/time'));
%time_spec=sorttime(time_spec);

%find temporal index
it=find(time_spec <= t_spec,1,'last');
%read all m modes
spec=hdf5read_slice_new(filename,pathname,[it-1 0 n0],[1 mmodes 1]);
%normalise to maximum value
spec_rel=spec/max(spec);
%find poloidal peak of mode
mspec_max=find(spec==max(spec),1,'last');
mspec_max=m(mspec_max);

%plot result
figure;
subplot(2,1,1)
hold on;
plot(m,spec,'bx-')
xlabel('m')
ylab=strcat(' ',varlabel,'_{',num2str(n0),',m}');
strcat(typename,ylab)
ylabel(strcat(typename,ylab))
title(strcat('t=',num2str(time_spec(it)),' peaks at m=',num2str(mspec_max)))
%define axis
if deltam > 0
axis_old=axis;
axis([maxis_min-deltam maxis_max+deltam axis_old(3) axis_old(4)])
%plot boundaries of filter
line([maxis_min maxis_min],[axis_old(3) axis_old(4)],'LineStyle','--','Color','black')
line([maxis_max maxis_max],[axis_old(3) axis_old(4)],'LineStyle','--','Color','black')
end


subplot(2,1,2)
hold on;
plot(m,spec_rel,'bx-')
xlabel('m')
ylabel(strcat(typename, ylab,'/ Max ',ylab))
title(strcat('t=',num2str(time_spec(it)),' peaks at m=',num2str(mspec_max)))

%define axis
if deltam > 0
axis([maxis_min-deltam maxis_max+deltam 0 1.1])
end

%plot boundaries of filter
line([maxis_min maxis_min],[0 1.1],'LineStyle','--','Color','black')
line([maxis_max maxis_max],[0 1.1],'LineStyle','--','Color','black')

sortie='OK';

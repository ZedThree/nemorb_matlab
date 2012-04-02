%filename='/fuslwl/scratch/ssaar/chease/ogyropsi.h5';
%filename='/fuslwb/scratch/ssaar/mast22807_025_2.equ_h5';
%filename='/home/ssaar/matlab/chease/ogyropsi.h5';
c1=load_chease_eq(filename);

figure(1);
  pcolor(c1.R([1:end 1],:),c1.Z([1:end 1],:),c1.J([1:end 1],:));
shading interp;
hold on
 plot(c1.R([1:end 1],end),c1.Z([1:end 1],end),'*-')
hold off
axis equal

figure(2);
subplot(2,1,1);
  plot(c1.psinorm,c1.dpdpsi);
subplot(2,1,2);
  plot(c1.psinorm,c1.fdfdpsi);


figure(3);
subplot(2,1,1);
  plot(c1.psinorm,c1.q);
subplot(2,1,2);
  plot(c1.psinorm,c1.p);



figure(4);
 plot(c1.R',c1.Z','-')
axis equal

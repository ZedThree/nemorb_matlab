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


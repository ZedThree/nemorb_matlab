shift=1;
global qmin qmax qs0 deltam nfilt1;
qmin=min(q);
qmax=max(q);
is=find(s_prof <= s_spec,1,'last');
ds_prof=s_prof(2)-s_prof(1);
ws1=(s_spec-s_prof(is))/ds_prof;
ws0=1.0-ws1;
qs0=ws0*q(is)+ws1*q(is+1);
%disp('WARNING:CORRECTED qs0 DUE TO BUG')
%qs0=1.42


%
if tmax > time(end)
disp('tmax modified because it was too big')
tmax=time(end)
end

%
if tmin > time(end)
disp('tmin modified because it was too big')
tmin=time(1)
end

%plot_spec_hdf5(filename,1,0,nfilt1,tmax);
%plot_spec_hdf5(filename,2,0,nfilt1,tmax);
%plot_spec_hdf5(filename,3,0,nfilt1,tmax);
%plot_spec_hdf5(filename,1,1,nfilt1,tmax);
%plot_spec_hdf5(filename,2,1,nfilt1,tmax);
%plot_spec_hdf5(filename,3,1,nfilt1,tmax);

nloc=find_max_tor_hdf5(filename,1,0,tmax)
plot_spec_hdf5(filename,1,0,nloc,tmax);
%nloc=find_max_tor_hdf5(filename,2,0,tmax)
plot_spec_hdf5(filename,2,0,nfilt1,tmax);
%nloc=find_max_tor_hdf5(filename,3,0,tmax)
plot_spec_hdf5(filename,3,0,nloc,tmax);
nloc=find_max_tor_hdf5(filename,1,1,tmax)
plot_spec_hdf5(filename,1,1,nloc,tmax);
%nloc=find_max_tor_hdf5(filename,2,1,tmax)
plot_spec_hdf5(filename,2,1,nfilt1,tmax);
%nloc=find_max_tor_hdf5(filename,3,1,tmax)
plot_spec_hdf5(filename,3,1,nloc,tmax);


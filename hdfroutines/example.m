time_potrz=hdf5read(filename,'/data/var2d/potrz/time');
r=hdf5read(filename,'/data/var2d/potrz/coord1');
z=hdf5read(filename,'/data/var2d/potrz/coord2');
srz=hdf5read(filename,'/data/var2d/potrz/srz');
nr=length(r);
nz=length(z);
%read first and last slice
potrz(:,:,1)=hdf5read_slice(filename,'/data/var2d/potrz/data',[0              0 0],[1 nz nr]);
potrz(:,:,2)=hdf5read_slice(filename,'/data/var2d/potrz/data',[nsteps_potrz-1 0 0],[1 nz nr]);

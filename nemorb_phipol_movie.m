function nemorb_phipol_movie(sim, step_size, ind)
%---------------
%
%Description
%Make a movie of phi on poloidal plane
%---------------
%
%Dummy arguments
%
%sim       = array containing the simulations data, created with loader.m
%step_size = step size: sim.time(end)/step_size = number of frames
%ind       = array containing indexes of simulations for the output of this function
%--------------
%
%Output
%
%A movie of the evolution of phi on the poloidal plane
%--------------

if exist('ind')==0
ind=1;
end

for i=1:length(ind) %loop on simulations
k=ind(i);
num_frames = sim(k).time(end)/step_size;
for t = 0:step_size:sim(k).time(end)
nemorb_phipolrz_time_hdf5(sim(k), t, true, false, ind);
end
makemovie
end


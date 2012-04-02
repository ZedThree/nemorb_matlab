function nemorb_profiles(sim, species, ind)
%################################
%nemorb_transport_transport(sim, ind)
%################################
%---------------
%
%Description
%
% Plot profiles of temperature, density, velocity
%
%---------------
%
%Dummy arguments
%
% sim     = array containing the simulations data, created with loader.m
% species = cell array of species names - enclose in curly braces!
% ind     = array containing indexes of simulations for the output of this function
%---------------
%
%Output
%
% Graphs of temperature, density, velocity profiles
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

for jj=1:length(species)

% Arrays may need transposing
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



end % Species loop
end % Simulation loop

function nemorb_plot_vexb(sim, ind)
%################################
%nemorb_plot_vexb(sim,ind)
%################################
%---------------
%
%Description
%
%Plot vexbnavg and vexbnmax against time
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
%---------------

global timelab
if exist('ind') == 0
ind=1;
end

pwd_old=pwd;
for i=1:length(ind) %Loop on simulations
k = ind(i);
cd(sim(k).path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read data
sim(k).vexbnavg = hdf5read(sim(k).filename,['/data/var0d/main_ions/' ...
		    'vexbnavg']);
sim(k).vexbnmax = hdf5read(sim(k).filename,['/data/var0d/main_ions/' ...
		    'vexbnmax']);

figure;
plot(sim(k).time,sim(k).vexbnavg)
hold on;
plot(sim(k).time,sim(k).vexbnmax,'k')
title(strcat('V_{ExB} for sim = ', sim(k).name))
xlabel(timelab{1})
legend('Avg. V_{ExB}', 'Max. V_{ExB}')

end %Loop on simulations

cd(pwd_old);
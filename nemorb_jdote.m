function nemorb_jdote(sim, info, species, ind)

%################################
%nemorb_jdote(sim,info,species,ind)
%################################
%---------------
%
%Description
%
%Plot growth rate of toroidal mode vs. time (intended for single modes)
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
ind = 1;
end

pwd_old = pwd;
for i=1:length(ind) % Loop on simulations
k = ind(i);
cd(sim(k).path);

background_efield_size = nemorb_genericornot(sim(k),info(k),'background_efield_size');
if background_efield_size == 0.0
continue
end

% Make things a bit easier
jdote_bg_tot  = sim(k).(species).jdote_bg_tot;
jdote_bg_par  = sim(k).(species).jdote_bg_par;
jdote_bg_grb  = sim(k).(species).jdote_bg_grb;
jdote_bg_curv = sim(k).(species).jdote_bg_curv;
jdote_bg_grp  = sim(k).(species).jdote_bg_grp;
jdote_bg_exb  = sim(k).(species).jdote_bg_exb;

% JdotE sans pertubed field contribution
jdote_bg_sans = jdote_bg_par + jdote_bg_grb + jdote_bg_curv + jdote_bg_grp;

% $$$ figure;
% $$$ plot(sim(k).time, jdote_bg_tot);
% $$$ hold on;
% $$$ plot(sim(k).time, jdote_bg_exb, 'r:');
% $$$ plot(sim(k).time, jdote_bg_sans, 'k');
% $$$ xlabel(timelab{1})
% $$$ ylabel('Power')
% $$$ title(sprintf('Rate of energy transfer for %s',sim(k).name))
% $$$ legend('Total','Perturbed E field','Sans','Location','NorthWest')

figure;
plot(sim(k).time, jdote_bg_sans, 'm');
hold on;
plot(sim(k).time, jdote_bg_par);
plot(sim(k).time, jdote_bg_curv, 'r');
plot(sim(k).time, jdote_bg_grp, 'g');
plot(sim(k).time, jdote_bg_grb, 'k');
xlabel(timelab{1})
ylabel('Power')
title(sprintf('Rate of energy transfer for %s',sim(k).name))
legend('Sans','Parallel','Curvature','Grad p','Grad B','Location','NorthWest')

figure;
plot(sim(k).time, (jdote_bg_curv +  jdote_bg_grb), 'm');
hold on;
plot(sim(k).time, jdote_bg_grb, 'k');
plot(sim(k).time, jdote_bg_curv, 'r');
xlabel(timelab{1})
ylabel('Power')
title(sprintf('Rate of energy transfer for %s',sim(k).name))
legend('Curvature','Grad B','Curv + grad B','Location','NorthWest')
% $$$ ratio1 = jdote_bg_curv./(jdote_bg_curv + jdote_bg_grb);
% $$$ ratio2 = jdote_bg_grb./(jdote_bg_curv + jdote_bg_grb);
% $$$ figure;
% $$$ plot(sim(k).time, ratio1);
% $$$ figure;
% $$$ plot(sim(k).time, ratio2,'r');

end


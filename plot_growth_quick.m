function plot_growth_quick(name)
% Stupid little script to just save growth rate graphs from the
% command line

% Open a simulation
[sim info]=nemorb_load({pwd},{name});
% Calculate the growth rates
nemorb_growth_convergance(sim,info,'deuterium',1,1);
% Save the graphs
saveas(1,'energy','fig');
saveas(2,'growth_rate','fig');
saveas(1,'energy','epsc2');
saveas(2,'growth_rate','epsc2');

exit

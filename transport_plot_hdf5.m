plot_ti=0;
plot_RoLTi=0;
plot_RoLNi=0;
plot_ni=0;
plot_etai=0;
plot_vpari=0;
plot_efluxi=0;
plot_chii=0;

plot_te=0;
plot_RoLTe=0;
plot_RoLNe=0;
plot_ne=0;
plot_etae=0;
plot_vpare=0;
plot_efluxe=0;
plot_chie=0;

if nsteps_1D > 1 %nsteps_1D gives a matlab error

%
%Plot evolution of temperature profile
%
if plot_ti==1
figure;
mytemp_i=double(mytemp_i);
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),mytemp_i(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of normalized ion temperature')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end
%

%
%Plot evolution of R/LT
%

if plot_RoLTi==1
RoLT=double(RoLT);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLT(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L_{T,i}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%
%Plot evolution of ni
%
if plot_ni==1
f_av_i=double(f_av_i);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),f_av_i(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of normalized ion density')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%
% Plot evolution of RoLNi
%
if plot_RoLNi == 1
RoLn=double(RoLn);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLn(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L{n,i}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%
% Plot evolution of etai
%

if plot_etai==1
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),etai(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \eta_i')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

if plot_vpari==1
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),v_par_av_hr(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of ion parallel momentum')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

%
%Plot evolution of Q
%
if plot_efluxi==1
figure;
efluxw_rad_i=double(efluxw_rad_i);
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),efluxw_rad_i(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of ion normalized heat flux')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

if plot_chii==1
chi_perp=double(chi_perp);
figure
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),chi_perp(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \chi_i')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

if nlelec > 0
%
%Plot evolution of temperature profile
%
if plot_te==1
figure;
mytemp_i_el=double(mytemp_i_el);
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),mytemp_i_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of normalized electron temperature')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end
%

%
%Plot evolution of R/LT
%

if plot_RoLTe==1
RoLT_el=double(RoLT_el);
RoLT_el_nop=double(RoLT_el_nop);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLT_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L_{T,e}')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLT_el_nop(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L_{T,e} without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

end

%
%Plot evolution of ne
%
if plot_ne==1
f_av_i_el_tot=double(f_av_i_el_tot);
f_av_i_el_tot_nop=double(f_av_i_el_tot_nop);

figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),f_av_i_el_tot(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of normalized electron density')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),f_av_i_el_tot_nop(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of normalized electron density without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

end

%
% Plot evolution of RoLNe
%
if plot_RoLNe == 1
RoLn_el=double(RoLn_el);
RoLn_el_nop=double(RoLn_el_nop);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLn_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L{n,e} with detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
%
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),RoLn_el_nop(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of R/L{n,e} without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

end

%
% Plot evolution of etae
%
if plot_etae==1
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),etae(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \eta_e')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

if plot_vpare==1
v_par_av_hr_el=double(v_par_av_hr_el);
v_par_av_hr_el_nop=double(v_par_av_hr_el_nop);
figure;
pcolor(time_1D,rad_prof1D_hr,v_par_av_hr_el)
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron parallel momentum with detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

figure;
pcolor(time_1D,rad_prof1D_hr,v_par_av_hr_el_nop)
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron parallel momentum without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])


end
%
%Plot evolution of Qe
%
if plot_efluxe==1
efluxw_rad_i_el=double(efluxw_rad_i_el);
efluxw_rad_i_el_nop=double(efluxw_rad_i_el_nop);
figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),efluxw_rad_i_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron normalized heat flux with detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

figure;
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),efluxw_rad_i_el_nop(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of electron normalized heat flux without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

end

if plot_chie==1
chi_perp_el=double(chi_perp_el);
chi_perp_el_nop=double(chi_perp_el_nop);
figure
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),chi_perp_el(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \chi_e with detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])

figure
pcolor(time_1D,rad_prof1D_hr(imin_1D:imax_1D),chi_perp_el_nop(imin_1D:imax_1D,:))
shading interp;
colorbar;
xlabel(timelab{1})
ylabel(radlab)
title('Evolution of \chi_e without detrapped')
set(gca,'Units','centimeters','Position',[ax ay aw ah])
end

end %nlelec
end %nsteps_1D

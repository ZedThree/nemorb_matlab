function [aspectratio tilt ellipse_strut] = nemorb_mode_aspectratio(sim, ...
						  show_fig, potsc, ind)
%



if exist('ind') == 0
ind=1;
end

for i=1:length(ind) % loop in simulations
k = ind(i);
cd(sim(k).path);

% get (R,Z) of potsc
r=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/rsc');
z=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/zsc');

% find global maximum and then local extrema at that radial surface
[row_max col_max]=find(potsc(:,:,k)==max(max(potsc(:,:,k))));
row_max=row_max(1);
[pks max_loc]=findpeaks(potsc(row_max,:,k));
[pks min_loc]=findpeaks(-potsc(row_max,:,k));

% if there's a minimum before a maximum, change sign of potential
% $$$ z_cutoff=min_loc(1);
% $$$ if (abs(max_loc(1)-min_loc(1)) > max_loc(1) || abs(max_loc(1)-min_loc(1)) ...
% $$$       == min_loc(1))
% $$$   potsc(:,:,k)=-potsc(:,:,k);
% $$$   z_cutoff=max_loc(2);
% $$$ end

extrema=sort([min_loc(1:5) max_loc(1:5)]);
z_cutoff=extrema(3);

% do interpolation
small_phi=griddata(r,z,potsc(:,:,k),r(1:end,1),z(row_max,1:z_cutoff),'cubic');

% take a small segment on the outboard midplane around first half-wavelength
%[row_min col_min]=find(potsc(:,:,k)==min(min(potsc(:,col_max:40,k))));
%[col_0]=find(potsc(row_max,col_max:col_min,k)<0,1);
s_start=1;
s_end=size(r,1);
chi_start=1;
%chi_end=col_max+col_0;
%chi_end=min_loc(find(min_loc > max_loc,1,'first'));
chi_end=z_cutoff;
%small_phi=potsc(s_start:s_end,chi_start:chi_end,k);
%small_r=r(s_start:s_end,chi_start:chi_end);
%small_z=z(s_start:s_end,chi_start:chi_end);
small_r=r(1:end,1);
small_z=z(row_max,1:z_cutoff);

% plot segment in background
loc_fig(k)=figure('Visible','off');
pcolor(small_r,small_z,small_phi);shading interp
%pcolor(ri,zi,potrz(:,:,k));shading interp
hold on
axis equal image

% get contour matrix of segment and half the max value of potsc
[c,h]=contour(small_r,small_z,small_phi);
if potsc(row_max,extrema(2),k) > 0
  midpoint=max(max(small_phi)).*0.5;
  [ii jj]=find(c>0 & c<midpoint,1,'last');
else
  midpoint=min(min(small_phi)).*0.5;
  [ii jj]=find(c<0 & c>midpoint,1,'last');
end
% get indices of first contour at midpoint level
kk=c(2,jj);
% contour at midpoint level for first mode half-wavelength
mode_half=c(:,jj+1:jj+kk);
% overplot selected contour, for debugging
plot(mode_half(1,:),mode_half(2,:),'+k','MarkerSize',8);

% get ellipse coefficients, overplot
ellipse_strut(k)=fit_ellipse(mode_half(1,:),mode_half(2,:),loc_fig(k));
title(sim(k).name);

aspectratio(k)=ellipse_strut(k).short_axis/ellipse_strut(k).long_axis;
tilt(k)=ellipse_strut(k).phi;

end


if show_fig==1
  set(loc_fig,'Visible','on');
else
  for i=1:length(ind)
    k=ind(i);
    close(loc_fig(k));
  end
end


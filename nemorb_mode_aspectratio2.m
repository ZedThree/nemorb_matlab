function [x_aspect delta_aspect ellipse_aspect ellipse_strut] = nemorb_mode_aspectratio2(sim, ...
						  show_fig, potsc, ind)
% [x_aspect delta_aspect ellipse_aspect ellipse_strut] = nemorb_mode_aspectratio2(sim, show_fig, potsc, ind)
% 
% Calculate aspect ratio of mode structure via three different methods:
% 	x_aspect     = aspect ratio of cross section through mode
% 		        maximum
%	delta_aspect = aspect ratio of total poloidal/radial
%			extents
%	ellipse_aspect = aspect ratio from fitted ellipse
%
% Inputs:
%	sim      = cell array of simulation data
%	show_fig = 1: display plots of potential and various fits
%	potsc    = poloidal cross section of potential
%	ind      = array of simulation indices

if exist('ind') == 0
ind=1;
end

pwd_old=pwd;

for i=1:length(ind) % loop in simulations
k = ind(i);
cd(sim(k).path);

% get (R,Z) of potsc
r=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/rsc');
z=hdf5read(sim(k).filename,'/data/var2d/generic/potsc/zsc');
nchi=size(r,2);

% find global maximum and then local extrema at that radial surface
[row_max col_max]=find(potsc(:,:,k)==max(max(potsc(:,:,k))));
row_max=row_max(1); % in case of multiple maxima
[max_pks max_loc]=findpeaks(potsc(row_max,:,k));
[min_pks min_loc]=findpeaks(-potsc(row_max,:,k));

% find extrema closest to outboard midplane
extrema=sort([min_loc max_loc]);
if (nchi - extrema(end)) < extrema(1)
  % mode is just below midplane
  small_z=[z(row_max,extrema(end-1):end) z(row_max,1:extrema(1))];
  midplane_mode=length(extrema);
else
  % mode is just above midplane
  small_z=[z(row_max,extrema(end):end) z(row_max,1:extrema(2))];
  midplane_mode=1;
end
r_extrema=r(row_max,extrema(midplane_mode));
z_extrema=z(row_max,extrema(midplane_mode));

% take a small segment on the outboard midplane around first half-wavelength
small_r=r(1:end,1);
try
  small_phi=griddata(r,z,potsc(:,:,k),small_r,small_z,'cubic');
catch err
  pottemp=potsc(:,:,k);
  F=TriScatteredInterp(r(:),z(:),pottemp(:));
  [qr qz]=meshgrid(small_r,small_z);
  small_phi=F(qr,qz);
end

% plot segment in background
if show_fig==1
  loc_fig(k)=figure('Visible','off');
  pcolor(small_r,small_z,small_phi);shading interp
  hold on
  axis equal image
end

% get contour matrix of segment at half the max value of potsc
midpoint=potsc(row_max,extrema(midplane_mode),k)*0.5;
[c num]=sensible_contour(small_r,small_z,small_phi,[midpoint midpoint]);

% there might multiple contours, so find one closest to our selected extrema
for ii=1:num
  avg_z(ii)=mean(c(ii).coords(2,:));
end
[dum ll]=min(abs(avg_z - z_extrema));
mode_contour=c(ll).coords;

% overplot selected contour
if show_fig==1
  plot(mode_contour(1,:),mode_contour(2,:),'+k','MarkerSize',8);
end

% aspect ratio of cross-section through extrema
rad_line=[small_r(1) small_r(end); z_extrema z_extrema];
z_line  =[r_extrema r_extrema; small_z(1) small_z(end)];
% find (R,Z) of lines through extrema
[x0 y0]=intersections(rad_line(1,:),rad_line(2,:),mode_contour(1,:),mode_contour(2,:));
[x1 y1]=intersections(z_line(1,:),z_line(2,:),mode_contour(1,:),mode_contour(2,:));

% error handling
if length(y1) < 2
  sprintf(['Not enough Z points to do cross-section aspect ratio for' ...
	   ' %s'],sim(k).name);
  return
end
if length(x0) < 2
  sprintf(['Not enough R points to do cross-section aspect ratio for' ...
	   ' %s'],sim(k).name);
  return
end
if length(y1) > 2
  y1=[max(y1) min(y1)];
end
if length(x0) > 2
  x0=[max(x0) min(x0)];
end
x_aspect(k)=abs(y1(2)-y1(1))/abs(x0(2)-x0(1));

% find delta-R,Z aspect ratio
deltaR(1)=min(mode_contour(1,:));
deltaR(2)=max(mode_contour(1,:));
deltaZ(1)=min(mode_contour(2,:));
deltaZ(2)=max(mode_contour(2,:));
delta_aspect(k)=abs(deltaZ(2)-deltaZ(1));%/abs(deltaR(2)-deltaR(1));

% delta_aspect = delta-Z through extrema
delta_aspect(k)=abs(y1(2)-y1(1));

% fit ellipse to mode structure
if show_fig==1
  ellipse_strut(k)=fit_ellipse(mode_contour(1,:),mode_contour(2,:),loc_fig(k));
else
  ellipse_strut(k)=fit_ellipse(mode_contour(1,:),mode_contour(2,:));
end
ellipse_aspect(k)=ellipse_strut(k).short_axis/ellipse_strut(k).long_axis;


% do plotting
if show_fig==1
  plot(x0, y0, 'm')
  plot(x1, y1, 'm')
  plot(r_extrema,z_extrema,'og')

  title(sim(k).name);
end

end % end simulation loop

% show graphs
if show_fig==1
  set(loc_fig,'Visible','on');
end

cd(pwd_old);
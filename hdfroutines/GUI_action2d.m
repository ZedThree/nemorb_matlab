%
% Callback functions for GUI_vis2d.m
%

function GUI_action2d(action)

global editHandle1 editHandle2 editHandle3 editHandle4 editHandle5...
       editHandle6 editHandle7 editHandle8 editHandle9 ...
       editHandle10 editHandle11 ...
       figureHandle1 figureHandle2 figureHandle3 figureHandle4 ...
       figureHandle5 figureHandle6 ...	 
       listboxHandle1 listboxHandle2 ...
       checkboxHandle1 checkboxHandle2 checkboxHandle3 ...
       checkboxHandle4 checkboxHandle5 ...
       popupmenuHandle1 popupmenuHandle2 ...
       popupmenuHandle3 popupmenuHandle4 ...
       axesHandle1 axesHandle2 axesHandle3 axesHandle4...
       axesHandle5 axesHandle6 ...
       sliderHandle1 ...
       buttonHandle1 buttonHandle2 buttonHandle3

global nx x_min x_max  xm dx ...
       ny y_min y_max  ym dy ...
       z_min z_max ...
       zx_min zx_max ...
       zy_min zy_max ...
       zxint_min zxint_max ...
       zyint_min zyint_max ...
       ind_x_min ind_x_max ind_y_min ind_y_max ...
       nx_skip ny_skip ...
       n_frames t_initial t_final tm i_frame time newY ...
       n_skip dt_pause ...
       field field_it ...
       plot_2dim plot_xcut plot_ycut plot_xint plot_yint ...
       plot_style shading_style ...
       x_text y_text z_text title_text...
       x_cut y_cut i_x_cut i_y_cut ...
       x_scale y_scale ...
       FileName VarName nc

switch action

case 'change_dir'

  dir_name = get(editHandle1, 'String');
  if ~isdir(dir_name),
    fprintf('!!! %s is not a directory !!!\n', dir_name)
    return; 
  end
  file_list = ls(dir_name);
  file_list = ['.. ', file_list];
  file_list = sort(str2cell(file_list));
  set(listboxHandle1, 'String', file_list, 'Value', [1])
   
case 'plot_2dim'

  plot_2dim = get(checkboxHandle1, 'Value');
  
  if plot_2dim == 1

    figureHandle2 = figure;
    set(figureHandle2, 'DeleteFcn', 'GUI_action2d close_2dim')
    axesHandle2 = axes('Position', [.1, .1, .8, .8], ...
    		       'XLim', [x_min x_max], ...
                       'YLim', [y_min y_max]);

    plot_field

  end
  
case 'read_plot_style'

  value = get(popupmenuHandle1, 'Value');
  plot_style = get(popupmenuHandle1, 'String');
  plot_style = plot_style{value};
  
  plot_field

case 'read_shading_style'

  value = get(popupmenuHandle2, 'Value');
  shading_style = get(popupmenuHandle2, 'String');
  shading_style = shading_style{value};

  plot_field

case 'read_XLim'

  eval(['x_lim =', get(editHandle8, 'String'), ';']);
  x_lim(1) = max(x_min, x_lim(1));
  x_lim(2) = min(x_max, x_lim(2));
  ind_x_min = interp1(xm, (1:nx), x_lim(1), 'nearest');
  ind_x_max = interp1(xm, (1:nx), x_lim(2), 'nearest');
  x_lim(1) = xm(ind_x_min);
  x_lim(2) = xm(ind_x_max);

  set(editHandle8, 'Enable', 'on', ...
    'String', ['[', num2str(x_lim(1), 4), ' ', num2str(x_lim(2), 4),']'])

  plot_field

case 'read_YLim'

  eval(['y_lim =', get(editHandle9, 'String'), ';']);
  y_lim(1) = max(y_min, y_lim(1));
  y_lim(2) = min(y_max, y_lim(2));
  ind_y_min = interp1(ym, (1:ny), y_lim(1), 'nearest');
  ind_y_max = interp1(ym, (1:ny), y_lim(2), 'nearest');
  y_lim(1) = ym(ind_y_min);
  y_lim(2) = ym(ind_y_max);

  set(editHandle9, 'Enable', 'on', ...
    'String', ['[', num2str(y_lim(1), 4), ' ', num2str(y_lim(2), 4),']'])

  plot_field
  
case 'read_xgrid_skip'

  nx_skip = get(editHandle10, 'String');
  nx_skip = str2num(nx_skip);
  nx_skip = abs(round(nx_skip));
  set(editHandle10, 'String', num2str(nx_skip));

  plot_field
  
case 'read_ygrid_skip'

  ny_skip = get(editHandle11, 'String');
  ny_skip = str2num(ny_skip);
  ny_skip = abs(round(ny_skip));
  set(editHandle11, 'String', num2str(ny_skip));

  plot_field
  
case 'close_2dim'

  plot_2dim = 0;
  set(checkboxHandle1, 'Value', [0])
  
case 'plot_xcut'

  plot_xcut = get(checkboxHandle2, 'Value');

  if plot_xcut == 1

    figureHandle3 = figure;
    set(figureHandle3, 'DeleteFcn', 'GUI_action2d close_xcut')
    axesHandle3 = axes('Position', [.1, .1, .8, .8], ...
                       'XLim', [y_min y_max]);

    x_half = (x_min + x_max)/2.;
    x_cut = x_half;
    set(editHandle2, 'Enable', 'on', 'String', num2str(x_half, 4))

    GUI_action2d read_xcut_position

  elseif plot_xcut == 0  
    GUI_action2d close_xcut
  end

case 'read_xcut_position'

  x_cut = get(editHandle2, 'String');
  x_cut = str2num(x_cut);
  x_cut = max(x_min, min(x_max, x_cut));
  i_x_cut = round(interp1(xm, (1:nx), x_cut));
  x_cut = xm(i_x_cut);

  set(editHandle2, 'String', num2str(x_cut, 4));

  zx_min = []; zx_max = []; 

  plot_field

case 'close_xcut'

  plot_xcut = 0;
  set(checkboxHandle2, 'Value', [0])
  set(editHandle2,     'Enable', 'off')
  zx_min = []; zx_max = []; 

case 'plot_ycut'

  plot_ycut = get(checkboxHandle3, 'Value');

  if plot_ycut == 1

    figureHandle4 = figure;
    set(figureHandle4, 'DeleteFcn', 'GUI_action2d close_ycut')
    axesHandle4 = axes('Position', [.1, .1, .8, .8], ...
                       'XLim', [x_min x_max]);

    y_half = (y_min + y_max)/2.;
    y_cut = y_half;
    set(editHandle3, 'Enable', 'on', 'String', num2str(y_half, 4))

    GUI_action2d read_ycut_position

  elseif plot_ycut == 0  
    GUI_action2d close_ycut
  end

case 'read_ycut_position'

  y_cut = get(editHandle3, 'String');
  y_cut = str2num(y_cut);
  y_cut = max(y_min, min(y_max, y_cut));
  i_y_cut = round(interp1(ym, (1:ny), y_cut));
  y_cut = ym(i_y_cut);

  set(editHandle3, 'String', num2str(y_cut, 4));

  zy_min = []; zy_max = [];

  plot_field

case 'close_ycut'

  plot_ycut = 0;
  set(checkboxHandle3, 'Value', [0])
  set(editHandle3,     'Enable', 'off')
  zy_min = []; zy_max = []; 

case 'plot_xint'

  plot_xint = get(checkboxHandle4, 'Value');

  if plot_xint == 1

    figureHandle5 = figure;
    set(figureHandle5, 'DeleteFcn', 'GUI_action2d close_xint')
    axesHandle5 = axes('Position', [.1, .1, .8, .8], ...
                       'XLim', [y_min y_max]);

    plot_field

  elseif plot_xint == 0
    GUI_action2d close_xint
  end

case 'close_xint'

  plot_xint = 0;
  set(checkboxHandle4, 'Value', [0])
  zxint_min = []; zxint_max = [];  

case 'plot_yint'

  plot_yint = get(checkboxHandle5, 'Value');

  if plot_yint == 1

    figureHandle6 = figure;
    set(figureHandle6, 'DeleteFcn', 'GUI_action2d close_yint')
    axesHandle6 = axes('Position', [.1, .1, .8, .8], ...
                       'XLim', [x_min x_max]);

    plot_field

  elseif plot_yint == 0
    GUI_action2d close_yint
  end

case 'close_yint'

  plot_yint = 0;
  set(checkboxHandle5, 'Value', [0])
  zyint_min = []; zyint_max = [];  

case 'read_x_scale'

  value = get(popupmenuHandle3, 'Value');
  x_scale = get(popupmenuHandle3, 'String');
  x_scale = x_scale{value};

  plot_field

case 'read_y_scale'

  value = get(popupmenuHandle4, 'Value');
  y_scale = get(popupmenuHandle4, 'String');
  y_scale = y_scale{value};

  plot_field

case 'run_forward'

  str = get(editHandle6, 'String');
  i_fr_initial = str2num(str);
  i_fr_initial = min(n_frames, max(1, i_fr_initial));

  set(buttonHandle3, 'UserData', [0])

  for i_frame = i_fr_initial:n_skip:n_frames,
	set(editHandle6, 'String', num2str(i_frame))
	GUI_action2d convert_frame

	value = get(buttonHandle3, 'UserData');

	if (value == 1) 
	  break
	end	

  end

case 'run_backward'

  str = get(editHandle6, 'String');
  i_fr_initial = str2num(str);
  i_fr_initial = min(n_frames, max(1, i_fr_initial));

  set(buttonHandle3, 'UserData', [0])

  for i_frame = i_fr_initial:-n_skip: 1,
	set(editHandle6, 'String', num2str(i_frame))
	GUI_action2d convert_frame

	value = get(buttonHandle3, 'UserData');

	if (value == 1) 
	  break
	end	

  end

case 'run_stop'

  set(buttonHandle3, 'UserData', [1])

case 'read_skip'

  n_skip = get(editHandle4, 'String');
  n_skip = str2num(n_skip);
  n_skip = abs(round(n_skip));
  set(editHandle4, 'String', num2str(n_skip));

case 'read_pause'

  dt_pause = get(editHandle5, 'String');
  dt_pause = str2num(dt_pause);
  dt_pause = abs(dt_pause);
  set(editHandle5, 'String', num2str(dt_pause));

case 'convert_frame'

  str = get(editHandle6, 'String');
  i_frame = str2num(str);
  i_frame = round(i_frame);
  i_frame = min(n_frames, max(1, i_frame));

  time = tm(i_frame);
  newY = (time - t_initial)/(t_final - t_initial); 

  plot_field

case 'convert_time'

  str = get(editHandle7, 'String');  
  time = str2num(str);
  time = min(t_final, max(t_initial, time));

  i_frame = round(interp1(tm, (1:n_frames), time));
  time = tm(i_frame);
  newY = (time - t_initial)/(t_final - t_initial);

  plot_field

case 'convert_slider'

  newY = get(sliderHandle1, 'Value');
  time = t_initial + newY*(t_final - t_initial);

  i_frame = round(interp1(tm, (1:n_frames), time));
  time = tm(i_frame);
  newY = (time - t_initial)/(t_final - t_initial);
  
  plot_field

case 'read_file'

  dir   = get(editHandle1, 'String');
  if strcmp(dir(end), '/'), dir = dir(1:end-1); end
  value = get(listboxHandle1, 'Value');
  FileName = get(listboxHandle1, 'String');
  FileName = deblank(FileName{value});
  
  if strcmp(FileName, '..');
    % Step back up in directory path
    is = strfind(dir, '/');
    if isempty(is) | is(end) == 1,
      FileName = '/';
    else
      FileName = dir(1:is(end)-1);
    end
  else
    FileName = [dir, '/', FileName];
    FileName = strrep(FileName, '//', '/');
  end
    
  if isdir(FileName),
    % If directory was chosen, step into directory
    set(editHandle1, 'String', FileName)
    GUI_action2d('change_dir')
    return
  end
  
  %-- Reading content of HDF5 file

  h5_info = hdf5info(FileName);

  if ~exist('h5_info'),
    fprintf('!!! File %s is not an HDF5 file!!!\n', FileName_1d)
    return
  end

  % find all 2-dim data (profiles)      

  cstr = vismd_GetVarNames(h5_info.GroupHierarchy, '2d');

  n_var2d = length(cstr);

  set(listboxHandle2, 'String', cstr)

case 'load_variable'

  set(checkboxHandle1,  'Enable', 'off', 'Value', [0])
  set(checkboxHandle2,  'Enable', 'off', 'Value', [0])
  set(checkboxHandle3,  'Enable', 'off', 'Value', [0])
  set(checkboxHandle4,  'Enable', 'off', 'Value', [0])
  set(checkboxHandle5,  'Enable', 'off', 'Value', [0])
  set(editHandle2,      'Enable', 'off')
  set(editHandle3,      'Enable', 'off')
  set(editHandle8,      'Enable', 'off')
  set(editHandle9,      'Enable', 'off')

  plot_2dim = 0;
  plot_xcut = 0;
  plot_ycut = 0;
  plot_xint = 0;
  plot_yint = 0;

  value = get(listboxHandle2, 'Value');
  VarName = get(listboxHandle2, 'String');
  VarName = VarName{value};

  % Save suffix of variable relative to species type
  % Suffix is of form '_ispecies'
  % This suffix is then used to load correct velocity mesh
  suffix = VarName(strfind(VarName, '_'):end);

  %--- load data from HDF5 file

  DataSetName = ['/data/var2d/', VarName, '/coord1'];
  [xm, attr] = hdf5read(FileName, DataSetName, 'ReadAttributes', true);
  x_text = attr.Value.Data;

  DataSetName = ['/data/var2d/', VarName, '/coord2'];
  [ym, attr] = hdf5read(FileName, DataSetName, 'ReadAttributes', true);
  y_text = attr.Value.Data;

  try
     tm = hdf5read(FileName, '/data/var2d/time');
  catch
     s  = lasterror;
     % If there isn't a time variable for all the 2d quantities,
     % there might be one in each directory.
     if findstr(s.message,'not an attribute or a dataset')==[]  
       rethrow(s);
     else 
        TimeName = ['/data/var2d/', VarName, '/time'];
        tm = hdf5read(FileName, TimeName);
     end
  end

  AttrName = ['/data/var2d/', VarName, '/title'];
  attr = hdf5read(FileName, AttrName);
  z_text = attr.Data;

  nx = length(xm);
  ny = length(ym);
  n_frames = length(tm);

  x_min = xm(1);
  x_max = xm(nx);
  dx    = xm(2) - xm(1);

  y_min = ym(1);
  y_max = ym(ny);
  dy    = ym(2) - ym(1);

  t_initial = tm(1);
  t_final   = tm(n_frames);

  z_min  = []; z_max  = [];
  zx_min = []; zx_max = [];
  zy_min = []; zy_max = [];
  zxint_min = []; zxint_max = [];
  zyint_min = []; zyint_max = [];

  i_frame = 1;
  time = t_initial;
  newY = 0.;

  set(checkboxHandle1,  'Enable', 'on')
  set(checkboxHandle2,  'Enable', 'on')
  set(checkboxHandle3,  'Enable', 'on')
  set(checkboxHandle4,  'Enable', 'on')
  set(checkboxHandle5,  'Enable', 'on')
  set(popupmenuHandle1, 'Enable', 'on')
  set(popupmenuHandle2, 'Enable', 'on')
  set(popupmenuHandle3, 'Enable', 'on')
  set(popupmenuHandle4, 'Enable', 'on')
  set(editHandle4,      'Enable', 'on')
  set(editHandle5,      'Enable', 'on')
  set(editHandle6,      'Enable', 'on')
  set(editHandle7,      'Enable', 'on')
  set(editHandle8,      'Enable', 'on')
  set(editHandle9,      'Enable', 'on')
  set(editHandle10,     'Enable', 'on')
  set(editHandle11,     'Enable', 'on')
  set(buttonHandle1,    'Enable', 'on')
  set(buttonHandle2,    'Enable', 'on')
  set(buttonHandle3,    'Enable', 'on')
  set(sliderHandle1,    'Enable', 'on', ...
	                'SliderStep', 1/(n_frames-1)*[1, 10])

  set(editHandle8, ...
      'String', ['[', num2str(x_min, 4), ' ', num2str(x_max, 4),']'])
  ind_x_min = 1;
  ind_x_max = nx;

  set(editHandle9, ...
      'String', ['[', num2str(y_min, 4), ' ', num2str(y_max, 4),']'])
  ind_y_min = 1;
  ind_y_max = ny;

  plot_field

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------

function plot_field

global editHandle1 editHandle2 editHandle3 editHandle4 editHandle5...
       editHandle6 editHandle7 editHandle8 editHandle9 ...
       figureHandle1 figureHandle2 figureHandle3 figureHandle4 ...
       figureHandle5 figureHandle6 ...	 
       listboxHandle1 listboxHandle2 ...
       checkboxHandle1 checkboxHandle2 checkboxHandle3 ...
       checkboxHandle4 checkboxHandle5 ...
       popupmenuHandle1 popupmenuHandle2 ...
       axesHandle1 axesHandle2 axesHandle3 axesHandle4...
       axesHandle5 axesHandle6 ...
       sliderHandle1 ...
       buttonHandle1 buttonHandle2 buttonHandle3

global nx x_min x_max  xm dx ...
       ny y_min y_max  ym dy ...
       z_min z_max ...
       zx_min zx_max ...
       zy_min zy_max ...
       zxint_min zxint_max ...
       zyint_min zyint_max ...
       nx_skip ny_skip ...
       ind_x_min ind_x_max ind_y_min ind_y_max ...
       n_frames t_initial t_final tm i_frame time newY ...
       n_skip dt_pause ...
       field field_it ...
       plot_2dim plot_xcut plot_ycut plot_xint plot_yint ...
       plot_style shading_style ...
       x_text y_text z_text title_text...
       x_cut y_cut i_x_cut i_y_cut ...
       x_scale y_scale ...
       FileName VarName nc

%--- reset slider and frame/time edit boxes

  set(editHandle6, 'String', num2str(i_frame))
  set(editHandle7, 'String', num2str(time, 4))
  set(sliderHandle1, 'Value', newY)

%--- find data of field relative to current time

%--- plot field
  VarName
  field_it = vismd_GetDataSet2(FileName, '2d', VarName, i_frame);
  field_it = field_it';

%--- title text

  is = strfind(FileName, '/');
  if length(is) <= 3,
    title_text = FileName;
  else
    title_text = ['...', FileName(is(end-2):end)];
  end
  title_text = ['File = '        , title_text, ...
	        '     time = '   , num2str(time, 4), ...
	        '     frame # = ',  num2str(i_frame)];

if plot_2dim == 1

  axes(axesHandle2)

  xp = xm(ind_x_min:nx_skip:ind_x_max);
  yp = ym(ind_y_min:ny_skip:ind_y_max);
  fp = double(field_it(ind_y_min:ny_skip:ind_y_max, ind_x_min:nx_skip:ind_x_max));
  
  eval([plot_style, '(xp, yp, fp)']);
  eval(['shading ', shading_style])

  z_min = min(z_min, min(min(fp)));
  z_max = max(z_max, max(max(fp)));
	
  if (z_min < z_max) 
    dlta_z = 0.2*(z_max-z_min);
    set(axesHandle2, 'ZLim', [z_min-dlta_z z_max+dlta_z])
  end

  h = colorbar('vert');
  hl = get(h, 'ylabel');
  set(hl, 'String', z_text, 'Rotation', [270]);

  xlabel(x_text); ylabel(y_text); zlabel(z_text)
  title(title_text, 'Interpreter', 'none')

  if plot_xcut == 1
	hold on
	plot([x_cut x_cut], [y_min y_max], 'k--', 'LineWidth', [4])
	hold off
  end

  if plot_ycut == 1
	hold on
	plot([x_min x_max], [y_cut y_cut], 'k--', 'LineWidth', [4])
	hold off
  end

end

if plot_xcut == 1
  
  switch x_scale,
  case 'lin',
    [zx_min, zx_max] = my_plot(ym, field_it(:, i_x_cut), ...
                               zx_min, zx_max, axesHandle3);
  case 'log', 
    my_semilogy(ym, field_it(:, i_x_cut), axesHandle3)
  end

  xlabel(y_text); ylabel(z_text)
  title(title_text, 'Interpreter', 'none')
  text(.1, .92, ['X-CUT = ', num2str(x_cut, 4)], ...
	'FontSize', [16], 'units', 'normalized')

end

if plot_ycut == 1
  
  switch y_scale,
  case 'lin',
    [zy_min, zy_max] = my_plot(xm, field_it(i_y_cut, :), ...
                            zy_min, zy_max, axesHandle4);
  case 'log', 
    my_semilogy(xm, field_it(i_y_cut, :), axesHandle4)
  end

  xlabel(x_text); ylabel(z_text)
  title(title_text, 'Interpreter', 'none')
  text(.1, .92, ['Y-CUT = ', num2str(y_cut, 4)], ...
	'FontSize', [16], 'units', 'normalized')

end

if plot_xint == 1
  
  % Trapezoidal integration:
  dist_xint = (sum(field_it(:, ind_x_min+1:ind_x_max-1), 2) ...
 	    + (field_it(:, ind_x_min) + field_it(:, ind_x_max))/2) ...
	    * dx;
%            / (ind_x_max - ind_x_min);


  switch x_scale,
  case 'lin',
    [zxint_min, zxint_max] = my_plot(ym, dist_xint, ...
                                     zxint_min, zxint_max, axesHandle5);
  case 'log', 
    my_semilogy(ym, dist_xint, axesHandle5)
  end

  xlabel(y_text); ylabel('\int dx')
  title(title_text, 'Interpreter', 'none')
  text(.1, .92, ['X-INT over x = [', ...
	num2str(xm(ind_x_min), 4), ', ', ...
	num2str(xm(ind_x_max), 4), ']'], ...
	'FontSize', [16], 'units', 'normalized')

end

if plot_yint == 1
  
  axes(axesHandle6)
	
  % Trapezoidal integration:
  dist_yint = (sum(field_it(ind_y_min+1:ind_y_max-1, :), 1) ...
	    + (field_it(ind_y_min, :) + field_it(ind_y_max, :))/2) ...
	    * dy;
%            / (ind_y_max - ind_y_min);

  switch y_scale,
  case 'lin',
    [zyint_min, zyint_max] = my_plot(xm, dist_yint, ...
                                     zyint_min, zyint_max, axesHandle6);
  case 'log', 
    my_semilogy(xm, dist_yint, axesHandle6)
  end

  xlabel(x_text); ylabel('\int dy')
  title(title_text, 'Interpreter', 'none')
  text(.1, .92, ['Y-INT over y = [', ...
	num2str(ym(ind_y_min), 4), ', ', ...
	num2str(ym(ind_y_max), 4), ']'], ...
	'FontSize', [16], 'units', 'normalized')

end

pause(dt_pause)

%----------------------------------------------------------------------
%----------------------------------------------------------------------

function cstr = str2cell(string)

%--- converts a string of words to a cell of strings for each word

iname = 0;
while ~isempty(string),
  iname = iname + 1;
  [name, string] = strtok(string);
  cstr(iname) = {name};
  string = deblank(string);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------

function my_semilogy(x, y, axesHandle)

axes(axesHandle)

y_plus = y;
y_plus(find(y_plus <= 0)) = NaN; 

y_minus = -y;
y_minus(find(y_minus <= 0)) =  NaN; 

% plotting positive values in blue, negative in red
semilogy(x, y_plus, 'b-')
hold on
semilogy(x, y_minus, 'r-')
hold off

%----------------------------------------------------------------------
%----------------------------------------------------------------------

function [y_min, y_max] = my_plot(x, y, y_min, y_max, axesHandle)

axes(axesHandle)

plot(x, y)

y_min = min([y_min, min(y)]);
y_max = max([y_max, max(y)]);

if (y_min < y_max),
  dlta_y = 0.2*(y_max-y_min);
  set(gca, 'YLim', [y_min-dlta_y y_max+dlta_y])
end

%
% Callback functions for GUI_vis1d.m
%

function GUI_action1d(action)

global editHandle1d_1 editHandle1d_2 editHandle1d_3 ...
       editHandle1d_4 editHandle1d_5 ...
       figureHandle1d_1 figureHandle1d_2 ...
       listboxHandle1d_1 listboxHandle1d_2 ...
       sliderHandle1d_1 ...
       buttonHandle1d_1 buttonHandle1d_2 buttonHandle1d_3 buttonHandle1d_4 ...
       popupmenuHandle1d_1 popupmenuHandle1d_2

global n_frames_1d t_initial_1d t_final_1d tm_1d ...
       i_frame_1d time_1d newY_1d ...
       n_skip_1d dt_pause_1d ...
       plot_obj n_hor n_ver plot_OK ...
       type_ylim y_scale ...
       FileName_1d ...
       title_filename_1d


switch action

case 'change_dir'

  dir_name = get(editHandle1d_1, 'String');
  if ~isdir(dir_name),
    fprintf('!!! %s is not a directory !!!\n', dir_name)
    return; 
  end
  file_list = ls(dir_name);
  file_list = ['.. ', file_list];
  file_list = sort(str2cell(file_list));
  set(listboxHandle1d_1, 'String', file_list, 'Value', [1])

case 'run_forward'

  str = get(editHandle1d_4, 'String');
  i_fr_initial = str2num(str) + 1;
  i_fr_initial = min(n_frames_1d, max(1, i_fr_initial));

  set(buttonHandle1d_3, 'UserData', [0])
	
  for i_frame_1d = i_fr_initial:n_skip_1d:n_frames_1d,
	set(editHandle1d_4, 'String', num2str(i_frame_1d))
	GUI_action1d convert_frame

	value = get(buttonHandle1d_3, 'UserData');

	if (value == 1) 
	  break
	end	

  end

case 'run_backward'

  str = get(editHandle1d_4, 'String');
  i_fr_initial = str2num(str) - 1;
  i_fr_initial = min(n_frames_1d, max(1, i_fr_initial));

  set(buttonHandle1d_3, 'UserData', [0])

  for i_frame_1d = i_fr_initial:-n_skip_1d: 1,
	set(editHandle1d_4, 'String', num2str(i_frame_1d))
	GUI_action1d convert_frame

	value = get(buttonHandle1d_3, 'UserData');

	if (value == 1) 
	  break
	end	

  end

case 'run_stop'

  set(buttonHandle1d_3, 'UserData', [1])

case 'read_skip'

  n_skip_1d = get(editHandle1d_2, 'String');
  n_skip_1d = str2num(n_skip_1d);
  n_skip_1d = abs(round(n_skip_1d));
  set(editHandle1d_2, 'String', num2str(n_skip_1d));

case 'read_pause'

  dt_pause_1d = get(editHandle1d_3, 'String');
  dt_pause_1d = str2num(dt_pause_1d);
  dt_pause_1d = abs(dt_pause_1d);
  set(editHandle1d_3, 'String', num2str(dt_pause_1d));

case 'convert_frame'
  %for the moment, use the first variable in the list...

  str = get(editHandle1d_4, 'String');
  i_frame_1d = str2num(str);
  i_frame_1d = round(i_frame_1d);
  i_frame_1d = min(length(plot_obj{1}.time), max(1, i_frame_1d));

  time_1d = plot_obj{1}.time(i_frame_1d);
  newY_1d = (time_1d - t_initial_1d)/(t_final_1d - t_initial_1d); 

  plot_profiles

case 'convert_time'

  str = get(editHandle1d_5, 'String');  
  time_1d = str2num(str);
  time_1d = min(t_final_1d, max(t_initial_1d, time_1d));
  newY_1d = (time_1d - t_initial_1d)/(t_final_1d - t_initial_1d);

  i_frame_1d = round(interp1(tm_1d, (1:n_frames_1d), time_1d));

  plot_profiles

case 'convert_slider'

  newY_1d = get(sliderHandle1d_1, 'Value');
  time_1d = t_initial_1d + newY_1d*(t_final_1d - t_initial_1d);

  i_frame_1d = round(interp1(tm_1d, (1:n_frames_1d), time_1d));

  plot_profiles

case 'read_y_limits'
  
  value = get(popupmenuHandle1d_1, 'Value');
  type_ylim = get(popupmenuHandle1d_1, 'String');
  type_ylim = type_ylim{value};  
  
  plot_profiles
  
case 'read_y_scale'

  value = get(popupmenuHandle1d_2, 'Value');
  y_scale = get(popupmenuHandle1d_2, 'String');
  y_scale = y_scale{value};

  plot_profiles
  
case 'read_file'
 
  % --- Read item chosen in listbox and identify whether it is a
  %     file or a directory 
  
  dir   = get(editHandle1d_1, 'String');
  if strcmp(dir(end), '/'), dir = dir(1:end-1); end
  value = get(listboxHandle1d_1, 'Value');
  FileName_1d= get(listboxHandle1d_1, 'String');
  FileName_1d= deblank(FileName_1d{value});
  
  if strcmp(FileName_1d, '..');
    % Step back up in directory path
    is = strfind(dir, '/');
    if isempty(is) | is(end) == 1,
      FileName_1d = '/';
    else
      FileName_1d = dir(1:is(end)-1);
    end
  else
    FileName_1d = [dir, '/', FileName_1d];
    FileName_1d = strrep(FileName_1d, '//', '/');
  end
    
  if isdir(FileName_1d),
    % If directory was chosen, step into directory
    set(editHandle1d_1, 'String', FileName_1d)
    GUI_action1d('change_dir')
    return
  end
  
  %--- File name text for title
  is = findstr(FileName_1d, '/');  
  if is <= 3,
    title_filename_1d = FileName_1d;
  else
    title_filename_1d = ['...', FileName_1d(is(end-2):end)];
  end


  %-- Reading content of HDF5 file
  h5_info = hdf5info(FileName_1d);

  if ~exist('h5_info'),
    fprintf('!!! File %s is not an HDF5 file!!!\n', FileName_1d)
    return
  end

  % find all 1-dim data (profiles)      
  cstr = vismd_GetVarNames(h5_info.GroupHierarchy, '1d');

  set(listboxHandle1d_2, 'String', cstr)
 
case 'plot_data'

  value = get(listboxHandle1d_2, 'Value');
  cstr = get(listboxHandle1d_2, 'String');
  cstr = cstr(value);

  if length(cstr)==0,
     fprintf('!!! No 1D data in this HDF5 file')
     return
  end
  
  n_var1d = length(cstr);
  plot_obj = cell(size(cstr));

  default_time_present = 1;
  time_found = 0;
  try
     [tm_1d, t_attr] = hdf5read(FileName_1d, '/data/var1d/time', 'ReadAttributes', true);
     t_text_1d = t_attr.Value.Data;
     time_found=1;
  catch
     fprintf('no overall time coord..');
     default_time_present = 0;
     s  = lasterror
     if findstr(s.message,'not an attribute or a dataset')==[]  
       rethrow(s);
     end
  end
 
  tm_1d_multi = {};

  for iv = 1:length(cstr),
    VarName = char(cstr(iv));

    AttrName = ['data/var1d/', VarName, '/title'];
    attr = hdf5read(FileName_1d, AttrName);
    data_text = attr.Data;

    DataSetName = ['/data/var1d/', VarName, '/coord1'];
    [xm, attr] = hdf5read(FileName_1d, DataSetName, 'ReadAttributes', true);
    x_text = attr.Value.Data;

    PlotOrder(iv) = hdf5read(FileName_1d, ['data/var1d/', VarName, '/PlotOrder']);

    % Try and read the individual time values..
    try
        TimeName = ['/data/var1d/', VarName, '/time'];
        [tm_1d_m, t_attr] = hdf5read(FileName_1d, TimeName,'ReadAttributes', true);      
        if (time_found==0),
	    time_found=1;
            tm_1d = tm_1d_m;
        end
    catch
        if (default_time_present==0),
	   fprintf('WARNING: no time variable available for a 1D quantity')
        end
	tm_1d_m = tm_1d;
        s  = lasterror;
        if findstr(s.message,'not an attribute or a dataset')==[]  
          rethrow(s);
        end
    end

    i_frame = 1;

    plot_obj{iv} = struct('name'     , VarName  , ...
                          'i_frame'  , i_frame  , ...
    	                  'data_text', data_text, ...
                          'time'     , tm_1d_m  , ...
			  'xm'       , xm       , ...
			  'x_text'   , x_text   , ...
    		          'ymin'     , []       , ...     
    		          'ymax'     , []            );
    tstart(iv) = tm_1d_m(1);
    tend(iv)   = tm_1d_m(end);
  end  

  t_text_1d = t_attr.Value.Data;
  
  %--- End reading HDF5 file

  % Order variables for plotting
  [dummy, ind_sort] = sort(PlotOrder);

  plot_obj = plot_obj(ind_sort);
  cstr = cstr(ind_sort);

  %--- Initializing

  % Frame id now not really used..
  % kept for compatibility.
  i_frame_id = 1;
  t_initial_1d = min(tstart);
  t_final_1d   = max(tend);

  n_hor = round(sqrt(n_var1d));  % number of horizontal plots
  n_ver = ceil(n_var1d/n_hor) ;  % number of vertical plots

  time_1d = t_initial_1d;
  newY_1d = 0.;

  figureHandle1d_2 = figure;
  set(figureHandle1d_2, 'DeleteFcn', 'GUI_action1d close_figure')

  plot_OK = 1;

  plot_profiles

  set(editHandle1d_2,   'Enable', 'on')
  set(editHandle1d_3,   'Enable', 'on')
  set(editHandle1d_4,   'Enable', 'on')
  set(editHandle1d_5,   'Enable', 'on')
  set(buttonHandle1d_1, 'Enable', 'on')
  set(buttonHandle1d_2, 'Enable', 'on')
  set(buttonHandle1d_3, 'Enable', 'on')
  %sensible when all have same time scale, otherwise problematic..
  n_frames_1d = length(plot_obj{1}.time);
  set(sliderHandle1d_1, 'Enable', 'on', ...
	                'SliderStep', 1/(n_frames_1d-1)*[1, 10])
  set(popupmenuHandle1d_1, 'Enable', 'on')
  set(popupmenuHandle1d_2, 'Enable', 'on')

 
case 'close_figure'

 plot_OK = 0;

end

%----------------------------------------------------------------------
%----------------------------------------------------------------------

function plot_profiles

global editHandle1d_1 editHandle1d_2 editHandle1d_3...
       editHandle1d_4 editHandle1d_5 ...
       figureHandle1d_1 figureHandle1d_2 ...
       listboxHandle1d_1 listboxHandle1d_2 ...
       sliderHandle1d_1 ...
       buttonHandle1d_1 buttonHandle1d_2 buttonHandle1d_3

global frames_1d t_initial_1d t_final_1d tm_1d ...
       i_frame_1d time_1d newY_1d ...
       n_skip_1d dt_pause_1d ...
       plot_obj n_hor n_ver plot_OK ...
       type_ylim y_scale ...
       FileName_1d ...
       title_filename_1d


%--- reset slider and frame/time edit boxes

  set(editHandle1d_4, 'String', num2str(i_frame_1d))
  set(editHandle1d_5, 'String', num2str(time_1d, 4))
  set(sliderHandle1d_1, 'Value', newY_1d)

%--- plot profiles
  if plot_OK == 1

    figure(figureHandle1d_2)

    n_var1d = length(plot_obj);
    for iv = 1:n_var1d
 
       i_frame_m = round(interp1(transpose(plot_obj{iv}.time), (1:length(plot_obj{iv}.time)), time_1d));

       xm = plot_obj{iv}.xm;
       y = vismd_GetDataSet2(FileName_1d, '1d', plot_obj{iv}.name, i_frame_m);

       axesHandle = subplot(n_ver, n_hor, iv);
       
       switch y_scale
        case 'lin',
         [plot_obj{iv}.ymin, plot_obj{iv}.ymax] = ...
             my_plot(xm, y, plot_obj{iv}.ymin, plot_obj{iv}.ymax, ...
                     axesHandle, type_ylim);
        case 'log',
         my_semilogy(xm, y, axesHandle)
       end

      xlabel(plot_obj{iv}.x_text); ylabel(plot_obj{iv}.data_text);
	
    end


    subplot(n_ver, n_hor, round(n_hor/2))
    title_text = ['File = '        , title_filename_1d   , ...
	          '     time = '   , num2str(time_1d, 4) , ...
	          '     frame # = ',  num2str(i_frame_1d)];
    title(title_text, 'Interpreter', 'none')

    pause(dt_pause_1d)

  end
 
%----------------------------------------------------------------------
%----------------------------------------------------------------------

function [y_min, y_max] = my_plot(x, y, y_min, y_max, axesHandle, type_ylim)

axes(axesHandle)

plot(x, y)

y_min = min([y_min, min(y)]);
y_max = max([y_max, max(y)]);

if strcmp(type_ylim, 'fixed') & (y_min < y_max),
  dlta_y = 0.2*(y_max-y_min);
  set(axesHandle, 'YLim', [y_min-dlta_y y_max+dlta_y])
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


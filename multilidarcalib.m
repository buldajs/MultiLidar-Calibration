function varargout = multilidarcalib(varargin)
% MULTILIDARCALIB MATLAB code for multilidarcalib.fig
%      MULTILIDARCALIB, by itself, creates a new MULTILIDARCALIB or raises the existing
%      singleton*.
%
%      H = MULTILIDARCALIB returns the handle to a new MULTILIDARCALIB or the handle to
%      the existing singleton*.
%
%      MULTILIDARCALIB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTILIDARCALIB.M with the given input arguments.
%
%      MULTILIDARCALIB('Property','Value',...) creates a new MULTILIDARCALIB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multilidarcalib_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multilidarcalib_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multilidarcalib

% Last Modified by GUIDE v2.5 24-Sep-2016 19:55:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multilidarcalib_OpeningFcn, ...
                   'gui_OutputFcn',  @multilidarcalib_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before multilidarcalib is made visible.
function multilidarcalib_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multilidarcalib (see VARARGIN)

% Choose default command line output for multilidarcalib
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes multilidarcalib wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multilidarcalib_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in about_bn.
function about_bn_Callback(hObject, eventdata, handles)
% hObject    handle to about_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
about_dlg;


% --- Executes on button press in exit_bn.
function exit_bn_Callback(hObject, eventdata, handles)
% hObject    handle to exit_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --- Executes on button press in search_scans_bn1.
function search_scans_bn1_Callback(hObject, eventdata, handles)
% hObject    handle to search_scans_bn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Query base name of lidar1 scans
temp=input('Enter base name of lidar1 scans (without numbers or suffix): ','s');
if(isempty(temp))
    disp('Invalid entry. Try again.');
    return;
end
handles.lidar1scanBaseName = temp;
% Query suffix name of lidar1 scans
temp=input('Enter suffix of lidar1 scan files ([] = "xyz"): ','s');
if(isempty(temp))
    temp = 'xyz';
end 
handles.lidar1scanSuffix=temp;
handles.lidar1_ind_active_range_images=[];

% Find lidar1 scan files in current directory whose indices match
% ind_active
d = dir([handles.lidar1scanBaseName '*.' handles.lidar1scanSuffix]);
nFiles =size(d,1);
if(~nFiles)
    disp(['Found no files named ' handles.lidar1scanBaseName ...
        '*.' handles.lidar1scanSuffix '! Check your input and try again. ']);
    return;
end

for i=1:nFiles
    % Strip path free version lidar1scanBaseName 
    s=findstr('/',handles.lidar1scanBaseName);
    if(~isempty(s))
        lidar1scanBaseName_nopath=handles.lidar1scanBaseName(max(s)+1:end);
    else
        lidar1scanBaseName_nopath=handles.lidar1scanBaseName;
    end
    s=findstr('\',lidar1scanBaseName_nopath);
    if(~isempty(s))
        lidar1scanBaseName_nopath=lidar1scanBaseName_nopath(max(s)+1:end);
    end
    % Extract number from filename
    s=findstr(d(i).name,lidar1scanBaseName_nopath);
    ss=findstr(d(i).name,['.' handles.lidar1scanSuffix]);
    filenum=str2num(d(i).name(s+length(lidar1scanBaseName_nopath):(ss-1)));
    % Add it to current active list
    handles.lidar1_ind_active_range_images=[handles.lidar1_ind_active_range_images filenum];
end
handles.lidar1_ind_active_range_images = sort(handles.lidar1_ind_active_range_images);
handles.lidar1_user_selected_planes = cell(1,length(handles.lidar1_ind_active_range_images));
disp('--- Loaded list of lidar1 active scan range image ---');
disp(int2str(handles.lidar1_ind_active_range_images));

% Update handles structure
guidata(hObject, handles);

set(handles.search_scans_bn2, 'Enable' , 'On');


% --- Executes on button press in search_scans_bn2.
function search_scans_bn2_Callback(hObject, eventdata, handles)
% hObject    handle to search_scans_bn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Query base name of scans
temp=input('Enter base name of lidar2 scans (without numbers or suffix): ','s');
if(isempty(temp))
    disp('Invalid entry. Try again.');
    return;
end

handles.lidar2scanBaseName = temp;
% Query suffix name of lidar2 scans
temp=input('Enter suffix of scan files ([] = "xyz"): ','s');
if(isempty(temp))
    temp = 'xyz';
end 
handles.lidar2scanSuffix=temp;
handles.ind_active_pairs=[];

% Find lidar2 scan files in current directory whose indices match
% ind_active
d=dir([handles.lidar2scanBaseName '*.' handles.lidar2scanSuffix]);
nFiles=size(d,1);
if(~nFiles)
    disp(['Found no files named ' handles.lidar2scanBaseName ...
        '*.' handles.lidar2scanSuffix '! Check your input and try again.']);
    return;
end

for i=1:nFiles
    % Strip path free version lidar2scanBaseName 
    s=findstr('/',handles.lidar2scanBaseName);
    if(~isempty(s))
        lidar2scanBaseName_nopath=handles.lidar2scanBaseName(max(s)+1:end);
    else
        lidar2scanBaseName_nopath=handles.lidar2scanBaseName;
    end
    s=findstr('\',lidar2scanBaseName_nopath);
    if(~isempty(s))
        lidar2scanBaseName_nopath=lidar2scanBaseName_nopath(max(s)+1:end);
    end
    % Extract number from filename
    s=findstr(d(i).name,lidar2scanBaseName_nopath);
    ss=findstr(d(i).name,['.' handles.lidar2scanSuffix]);
    filenum=str2num(d(i).name(s+length(lidar2scanBaseName_nopath):(ss-1)));
    % If number is in lidar1_ind_active_range_images from lidar1 selected
    % result
    if(~isempty(find(handles.lidar1_ind_active_range_images==filenum)))
        % Add it to current active list
        handles.ind_active_pairs=[handles.ind_active_pairs filenum];
    end
end

handles.ind_active_pairs=sort(handles.ind_active_pairs);
handles.lidar2_user_selected_planes=cell(1,length(handles.ind_active_pairs));

%update lidar1_ind_active_range_images and lidar1_user_selected_planes
handles.lidar1_ind_active_range_images = sort(handles.ind_active_pairs);
handles.lidar1_user_selected_planes = cell(1,length(handles.ind_active_pairs));

disp('--- Loaded list of active lidar1 scan/lidar2 scan pairs ---');
disp(int2str(handles.ind_active_pairs));


% Update handles structure
guidata(hObject, handles);

set(handles.select_lidar1_planes_bn, 'Enable','On');
set(handles.select_lidar2_planes_bn, 'Enable','On');

% --- Executes on button press in select_lidar1_planes_bn.
function select_lidar1_planes_bn_Callback(hObject, eventdata, handles)
% hObject    handle to select_lidar1_planes_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Preview a set of lidar range images in separate window
ind_preview=input('Enter indices of lidar1 range images to preview ([]=all active): ');
if(isempty(ind_preview))
    if(isempty(handles.lidar1_ind_active_range_images))
        disp('No active images to preview');
        return;
    end
    ind_preview=handles.lidar1_ind_active_range_images;
end

% Create cell array containing valid filenames for preview
preview_filenames={};
nvalid=0;
for i=1:length(ind_preview)
    % Check to see if it is a valid lidar1 range image file
    d=dir([handles.lidar1scanBaseName int2str(ind_preview(i)) '.' handles.lidar1scanSuffix]);
    if(~size(d,1))
        continue;
    end
    nvalid=nvalid+1;
    preview_filenames{nvalid}=[handles.lidar1scanBaseName int2str(ind_preview(i)) '.' handles.lidar1scanSuffix];
end

% Construct a cell containing info about previously lidar1 selected planes
old_lidar1_selected_planes=cell(1,length(ind_preview));
for i=1:length(ind_preview)
    arrayPos=find(handles.lidar1_ind_active_range_images==ind_preview(i));
    old_lidar1_selected_planes{i}=...
        handles.lidar1_user_selected_planes{arrayPos};
end

disp('Starting lidar1 range image viewer... please wait');
% View range images in "select planes" mode
[h, temp]=view_range_image(preview_filenames,1,old_lidar1_selected_planes,1,10);
delete(h);
clear old_lidar1_selected_planes;

% If user pushed the "Cancel" button
if(isempty(temp))
    return;
end

% Update the contents of handles.lidar1_user_selected_plane with the new info
for i=1:length(temp)
    if(~isempty(temp{i}))
        arrayPos=find(handles.lidar1_ind_active_range_images==ind_preview(i));
        handles.lidar1_user_selected_planes{arrayPos}=temp{i};
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in select_lidar2_planes_bn.
function select_lidar2_planes_bn_Callback(hObject, eventdata, handles)
% hObject    handle to select_lidar2_planes_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Preview a set of range images in separate window
ind_preview=input('Enter indices of lidar2 range images to preview ([]=all active): ');
if(isempty(ind_preview))
    if(isempty(handles.ind_active_pairs))
        disp('No active images to preview');
        return;
    end
    ind_preview=handles.ind_active_pairs;
end

% Create cell array containing valid filenames for preview
preview_filenames={};
nvalid=0;
for i=1:length(ind_preview)
    % Check to see if it is a valid range image file
    d=dir([handles.lidar2scanBaseName int2str(ind_preview(i)) '.' handles.lidar2scanSuffix]);
    if(~size(d,1))
        continue;
    end
    nvalid=nvalid+1;
    preview_filenames{nvalid}=[handles.lidar2scanBaseName int2str(ind_preview(i)) '.' handles.lidar2scanSuffix];
end

% Construct a cell containing info about previously selected planes
lidar2_old_selected_planes=cell(1,length(ind_preview));
for i=1:length(ind_preview)
    arrayPos=find(handles.lidar1_ind_active_range_images==ind_preview(i));
    lidar2_old_selected_planes{i}=...
        handles.lidar2_user_selected_planes{arrayPos};
end


disp('Starting lidar2 range image viewer... please wait');
% View range images in "select planes" mode
[h, temp]=view_range_image(preview_filenames,1,lidar2_old_selected_planes,1,10);
delete(h);
clear lidar2_old_selected_planes;

% If user pushed the "Cancel" button
if(isempty(temp))
    return;
end
% Update the contents of handles.user_selected_plane with the new info
for i=1:length(temp)
    if(~isempty(temp{i}))
        arrayPos=find(handles.lidar1_ind_active_range_images==ind_preview(i));
        handles.lidar2_user_selected_planes{arrayPos}=temp{i};
    end
end
% Update handles structure
guidata(hObject, handles);


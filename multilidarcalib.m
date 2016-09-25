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

% Last Modified by GUIDE v2.5 25-Sep-2016 14:41:05

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

set(handles.select_lidar2_planes_bn, 'Enable','On');


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

set(handles.calibrate_bn, 'Enable','On');


% --- Executes on button press in calibrate_bn.
function calibrate_bn_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lidar1 part:
% For each laser scan corresponding to the used images
% Select the 3d points corresponding to the target from the range image
% Fit a plane to the points and get [theta,alpha] wrt the laser 
lidar1_n_scans=length(handles.lidar1_ind_active_range_images);

lidar1_thetal=[]; lidar1_alphal=[];
lidar1_planePoints=[];
lidar1_n_planePoints=[];
%weight=[];

for j=1:lidar1_n_scans
    i=find(handles.lidar1_ind_active_range_images==handles.ind_active_pairs(j));
    fprintf(1,'Processing %s\n',[handles.lidar1scanBaseName ...
        int2str(handles.lidar1_ind_active_range_images(j)) '.' handles.lidar1scanSuffix]);
    % theta and alpha in lidar1 frame
    lidar1_thetal=[lidar1_thetal handles.lidar1_user_selected_planes{i}.theta];
    lidar1_alphal=[lidar1_alphal handles.lidar1_user_selected_planes{i}.alpha];
    % inliers in lidar1 selected plane points
    lidar1_planePoints=[lidar1_planePoints; handles.lidar1_user_selected_planes{i}.inliers];
    % number of inliers in selected plane points
    lidar1_n=size(handles.lidar1_user_selected_planes{i}.inliers,1);
    lidar1_n_planePoints=[lidar1_n_planePoints ; lidar1_n];
    % Each point is weighted inversely to the number of inliers selected 
    % in each selected plane
    %weight=[weight; ones(n,1)./n];
    fprintf(1,'Median error was %f\n',handles.lidar1_user_selected_planes{i}.e);
end


% Lidar2 part:
% For each laser scan corresponding to the used images
% Select the 3d points corresponding to the target from the range image
% Fit a plane to the points and get [theta,alpha] wrt the laser 
lidar2_n_scans=length(handles.ind_active_pairs);

lidar2_thetal=[]; lidar2_alphal=[];
lidar2_planePoints=[];
lidar2_n_planePoints=[];
%weight=[];

for j=1:lidar2_n_scans
    i=find(handles.lidar1_ind_active_range_images==handles.ind_active_pairs(j));
    fprintf(1,'Processing %s\n',[handles.lidar2scanBaseName ...
        int2str(handles.ind_active_pairs(j)) '.' handles.lidar2scanSuffix]);
    % theta and alpha in lidar2 frame
    lidar2_thetal=[lidar2_thetal handles.lidar2_user_selected_planes{i}.theta];
    lidar2_alphal=[lidar2_alphal handles.lidar2_user_selected_planes{i}.alpha];
    % inliers in lidar2 selected plane points
    lidar2_planePoints=[lidar2_planePoints; handles.lidar2_user_selected_planes{i}.inliers];
    % number of inliers in selected plane points
    lidar2_n=size(handles.lidar2_user_selected_planes{i}.inliers,1);
    lidar2_n_planePoints=[lidar2_n_planePoints ; lidar2_n];
    % Each point is weighted inversely to the number of inliers selected 
    % in each selected plane
    %weight=[weight; ones(n,1)./n];
    fprintf(1,'Median error was %f\n',handles.lidar2_user_selected_planes{i}.e);
end

disp('-- Optimization: Stage I');
% Computing the best transformation:
% We proceed by performing the minimization in a iterated two-stage 
% process.
% In the first stage, we find the translation that minimizes the
% difference in distance from the camera origin to each plane represented
% in the camera system and laser system.
% In the second stage, we find the best rotation that minimizes the angular
% difference between the normal from the origin to corresponding planes

% Assume planes params wrt lidar1 and lidar2 are [lidar1_thetal,lidar1_alphal] and 
% [lidar2_thetal,lidar2_alphal]. thetax is 3-by-n_planes , alphax is 1-by-n_planes
% 
t1=zeros(3,1); R1=eye(3);

% Closed form solution
t_est=inv(lidar1_thetal*lidar1_thetal')*lidar1_thetal*(lidar1_alphal-lidar2_alphal)';

% Adjust by t_estt_start=mean(thetal.*repmat(lidar2_alphal,3,1),2) ...

t1=t1+t_est;

fprintf(1,'Computed RMS error in distance to planes: %f\n',...
    computeRMSDiffDistanceToPlanes(t_est,lidar1_thetal,lidar1_alphal,lidar2_thetal,lidar2_alphal));

% Computing best rotation
% lidar1_normals(Q)= R * lidar2_normals(P)
[U,S,V]=svd(lidar2_thetal*lidar1_thetal');
R=V*U';
if(det(R)<0)
    fprintf(1,'*** WARNING : Found det(R) < 0 => optimal solution includes a reflection\n');
    resp = input('Allow reflection (y/n) [n]? ');    
    if ( ~isempty(resp) && lower(resp(1)) == 'y')
        R = V * diag([ ones(size(V,2)-1, 1) ; -1]) * U';
    else
        error('Exiting due to unexpect point configuration');
    end
end
R1=R;

disp('-- Optimization: Stage II');
%% Stage II of optimization
% Use the computed estimates of R1 and t1 in a second optimization
% stage where the sum of sqr distance of the polygon vertices to each plane
% is minimized
% Find quaternion corresponding to computed R1
q=rot2quat(R1);
par=[q ; t1];

cost=computeRMSWeightedDistVerticesToPlanes(par,lidar2_planePoints,lidar2_n_planePoints,...
    lidar1_thetal,lidar1_alphal);
fprintf(1,'Initial RMS distance of points to planes: %f\n',cost);
disp('... running non-linear optimization routine...');

% Uses optimization toolbox
% [par_est, fval, exitFlag,output]=fminunc(@computeDistVerticesToPlanes,...
%     par,[],planePoints,n_planePoints,thetac,alphac);
[par_est, fval, exitFlag,output]=fminsearch(...
    @computeRMSWeightedDistVerticesToPlanes,...
    par,[],lidar2_planePoints,lidar2_n_planePoints,lidar1_thetal,lidar1_alphal);

cost=computeRMSWeightedDistVerticesToPlanes(par_est,lidar2_planePoints,...
    lidar2_n_planePoints,lidar1_thetal,lidar1_alphal);
fprintf(1,'RMS distance of points to planes after search: %f\n',cost);

% Convert back to rotation matrix
R2=quat2rot(par_est(1:4));
t2=par_est(5:7);
%T=[quat2rot(par_est(1:4)) par_est(5:7)];

handles.R1=R1;
handles.t1=t1;
handles.R2=R2;
handles.t2=t2;

disp('-- Result: ');
[R1 t1]
[R2 t2]

% Take preventive measure in case user tries to quit without saving
handles.user_saved_results = 0;

% Update handles structure
guidata(hObject, handles);


% Enable save calib file  buttons
set(handles.save_calib_file_bn,'Enable','On');

return;

% =========================================================================
% Function takes as input a 7-vector [q t] corresp to quaterion
% representation of rotation, and translation t and computes the sum of
% square distances of the points, planeVertices, to the 3d planes.
% planeVertices are in the laser coordinate frame and contain 4*n_planes
% points (4 points per plane)
function cost = computeRMSWeightedDistVerticesToPlanes(par, ...
    planePoints,n_planePoints,thetac,alphac)

q = par(1:4);
t = par(5:7);
R = quat2rot(q);
cost=0; count=0;
cn_planePoints=cumsum(n_planePoints);
for i=1:length(n_planePoints)    
    %vert = planeVertices((4*i-3):(4*i),:);
    if(i==1)
        vert = planePoints(1:cn_planePoints(i),:);
    else
        vert = planePoints((cn_planePoints(i-1)+1):cn_planePoints(i),:);
    end
    theta = thetac(:,i);
    % Sum of square distances of points transformed to camera frame to
    % planes in camera frame
    cost=cost+mean((theta'*(R*vert'+repmat(t,1,n_planePoints(i)))-alphac(i)).^2);
end
cost=sqrt(cost./length(n_planePoints));
return;

%-------------------------------------------------------------------------    
function cost = computeRMSDiffDistanceToPlanes(pc,thetac,alphac,thetal,alphal)

%t_1= (abs(pc'*thetac - alphac) - alphal).^2 ;
%t_2= (abs(-pc'*thetal - alphal) - alphac).^2;
%t_3=max(t_1,t_2);
%cost=sum(t3(t3<median(t3)));
cost = sqrt(mean( (abs(pc'*thetac - alphac) - alphal).^2 ));
return;


% --- Executes on button press in save_calib_file_bn.
function save_calib_file_bn_Callback(hObject, eventdata, handles)
% hObject    handle to save_calib_file_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tag=input('Enter name tag for calibration file: ','s');
comment=input('Enter optional comment ([]=none): ','s');
% Save results of optim 1
fid=fopen([tag '_calib_1.m'],'wt');
fprintf(fid,'%% Lidar1 to Lidar2 calibration parameters (I optim stage) \n');
if (~isempty(comment))
    fprintf(fid,['%% File comment: ' comment '\n']);
end
fprintf(fid,['%% ' datestr(now) '\n']);
fprintf(fid,'%% \n');
fprintf(fid,'%% Transformation matrix specifies laser coordinate frame\n');
fprintf(fid,'%% in the reference frame of the camera\n');
fprintf(fid,'%% \n');
fprintf(fid,'%%-- Translation vector (t)\n');
fprintf(fid,'t = [ %5.6f ; %5.6f ; %5.6f ]\n',...
    handles.t1(1),handles.t1(2),handles.t1(3));
fprintf(fid,'%%-- Rotation matrix (R)\n');
fprintf(fid,'R = ...\n[ %5.6f  %5.6f  %5.6f ;...\n  %5.6f  %5.6f  %5.6f ;...\n  %5.6f  %5.6f  %5.6f ]\n',...
    handles.R1(1,1),handles.R1(1,2),handles.R1(1,3),...
    handles.R1(2,1),handles.R1(2,2),handles.R1(2,3),...
    handles.R1(3,1),handles.R1(3,2),handles.R1(3,3));
fclose(fid);
disp(['Saved ' tag '_calib_1.m']);

% Save results of optim 2
fid=fopen([tag '_calib_2.m'],'wt');
fprintf(fid,'%% Laser to Camera calibration parameters (II optim stages) \n');
if (~isempty(comment))
    fprintf(fid,['%% File comment: ' comment '\n']);
end
fprintf(fid,['%% ' datestr(now) '\n']);
fprintf(fid,'%% \n');
fprintf(fid,'%% Transformation matrix specifies laser coordinate frame\n');
fprintf(fid,'%% in the reference frame of the camera\n');
fprintf(fid,'%% \n');
fprintf(fid,'%%-- Translation vector (t)\n');
fprintf(fid,'t = [ %5.6f ; %5.6f ; %5.6f ]\n',...
    handles.t2(1),handles.t2(2),handles.t2(3));
fprintf(fid,'%%-- Rotation matrix (R)\n');
fprintf(fid,'R = ...\n[ %5.6f  %5.6f  %5.6f ;...\n  %5.6f  %5.6f  %5.6f ;...\n  %5.6f  %5.6f  %5.6f ]\n\n',...
    handles.R2(1,1),handles.R2(1,2),handles.R2(1,3),...
    handles.R2(2,1),handles.R2(2,2),handles.R2(2,3),...
    handles.R2(3,1),handles.R2(3,2),handles.R2(3,3));
fclose(fid);
disp(['Saved ' tag '_calib_2.m']);

handles.user_saved_results = 1;
% Update handles structure
guidata(hObject, handles);

function varargout = about_dlg(varargin)
% ABOUT_DLG MATLAB code for about_dlg.fig
%      ABOUT_DLG, by itself, creates a new ABOUT_DLG or raises the existing
%      singleton*.
%
%      H = ABOUT_DLG returns the handle to a new ABOUT_DLG or the handle to
%      the existing singleton*.
%
%      ABOUT_DLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUT_DLG.M with the given input arguments.
%
%      ABOUT_DLG('Property','Value',...) creates a new ABOUT_DLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before about_dlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to about_dlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help about_dlg

% Last Modified by GUIDE v2.5 24-Sep-2016 15:21:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @about_dlg_OpeningFcn, ...
                   'gui_OutputFcn',  @about_dlg_OutputFcn, ...
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


% --- Executes just before about_dlg is made visible.
function about_dlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to about_dlg (see VARARGIN)

% Choose default command line output for about_dlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes about_dlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = about_dlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ok_bn.
function ok_bn_Callback(hObject, eventdata, handles)
% hObject    handle to ok_bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

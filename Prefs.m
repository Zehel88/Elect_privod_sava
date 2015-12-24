function varargout = Prefs(varargin)
% PREFS MATLAB code for Prefs.fig
%      PREFS, by itself, creates a new PREFS or raises the existing
%      singleton*.
%
%      H = PREFS returns the handle to a new PREFS or the handle to
%      the existing singleton*.
%
%      PREFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREFS.M with the given input arguments.
%
%      PREFS('Property','Value',...) creates a new PREFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Prefs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Prefs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Prefs

% Last Modified by GUIDE v2.5 24-Dec-2015 18:43:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Prefs_OpeningFcn, ...
                   'gui_OutputFcn',  @Prefs_OutputFcn, ...
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


% --- Executes just before Prefs is made visible.
function Prefs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Prefs (see VARARGIN)

% Choose default command line output for Prefs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%% Начальные данные
if exist('Prefs.mat')==2
load('Prefs.mat');
% Загуржаем ранее использованные параметры
set(handles.uitable1,'DaTa',Prefs.params);
end


% UIWAIT makes Prefs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Prefs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Кнопка "применить"
function pushbutton1_Callback(hObject, eventdata, handles)
%% Применить изменения
% Сохраняем новые значения параметров
Prefs.params=get(handles.uitable1,'DaTa');

if isempty(Prefs.params{3,4})==0
save('Prefs.mat','Prefs');

close('Prefs');
else
    errordlg('Задайте параметры!!!');
end






% --- Кнопка "Отмена"
function pushbutton2_Callback(hObject, eventdata, handles)
%% Отменить изменения
Prefs.params=get(handles.uitable1,'DaTa');
if isempty(Prefs.params{3,4})==0
close('Prefs');
else
   errordlg('Задайте параметры!!!'); 
end

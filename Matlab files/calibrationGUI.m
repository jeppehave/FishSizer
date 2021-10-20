function varargout = calibrationGUI(varargin)
% CALIBRATIONGUI MATLAB code for calibrationGUI.fig
%      CALIBRATIONGUI, by itself, creates a new CALIBRATIONGUI or raises the existing
%      singleton*.
%
%      H = CALIBRATIONGUI returns the handle to a new CALIBRATIONGUI or the handle to
%      the existing singleton*.
%
%      CALIBRATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATIONGUI.M with the given input arguments.
%
%      CALIBRATIONGUI('Property','Value',...) creates a new CALIBRATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibrationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibrationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibrationGUI

% Last Modified by GUIDE v2.5 08-Mar-2021 14:08:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibrationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @calibrationGUI_OutputFcn, ...
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


% --- Executes just before calibrationGUI is made visible.
function calibrationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibrationGUI (see VARARGIN)

% Choose default command line output for calibrationGUI


global pic

currentFolder = getappdata(0,'folder1')

rstart=getappdata(0,'restart');
if rstart==1
    handles.calpath = getappdata(0,'cpath');
    oldfolder=cd(handles.calpath);
end

testpath=getappdata(0,'tpath')

if testpath~=0
    oldfolder=cd(testpath);
end

[picfile,picpath]= uigetfile( ...
{'*.jpg;*.jpeg;*.tif;*.tiff;*.bmp',...
    'Image files (*.jpg,*.jpeg,*.tif,*.tiff,*.bmp)';
   '*.*',  'All Files (*.*)'}, ...
   'Select calibration image file');

if picpath==0
    cd(currentFolder)
  return
else


pic=imread(fullfile(picpath,picfile));






if rstart==1 | testpath~=0
    cd(oldfolder)
end

axes(handles.axes1)
imshow(pic)

calllength=1;
set(handles.setlength,'string',sprintf('%f',calllength));

calcheck=0;
setappdata(0,'cc',calcheck);

setappdata(0,'tpath',picpath);

handles.calllength=calllength;

end
cd(currentFolder)
what
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes calibrationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calibrationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function setlength_Callback(hObject, eventdata, handles)
% hObject    handle to setlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setlength as text
%        str2double(get(hObject,'String')) returns contents of setlength as a double


% --- Executes during object creation, after setting all properties.
function setlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scale=handles.pmm;
calcheck=1;

oldfolder = getappdata(0,'folder1');

rstart=getappdata(0,'restart');
if rstart==1
    handles.calpath = getappdata(0,'cpath');
    oldfolder=cd(handles.calpath);
end

testpath=getappdata(0,'tpath')

if testpath~=0
    oldfolder=cd(testpath);
end

[calfile,calpath] = uiputfile('calibrationScale.mat','Save calibration')
save(fullfile(calpath,calfile),'scale')

oldfolder = getappdata(0,'folder1');
 cd(oldfolder)
 
if calpath==0
  % user pressed cancel
  return
end


if rstart==1 | testpath~=0
    cd(oldfolder)
end


setappdata(0,'pmm1',handles.pmm);
setappdata(0,'cc',calcheck);
setappdata(0,'cpath',calpath);
uiresume
guidata(hObject, handles);

cd(getappdata(0,'folder1'))

close calibrationGUI

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%pic=handles.pic
global pic;
axes(handles.axes1);
yl = ylim
xl = xlim
imshow(pic)
ylim(yl)
xlim(xl)
calllength=str2double(get(handles.setlength,'String'));
[x,y]=ginput(2);
dist1 = hypot(diff(x), diff(y));
pmm=calllength/dist1;
handles.pmm=pmm;

line(x,y,'LineWidth',2,'color','red')
guidata(hObject, handles);









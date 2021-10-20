function varargout = loadFishSizer(varargin)
% LOADFISHSIZER MATLAB code for loadFishSizer.fig
%      LOADFISHSIZER, by itself, creates a new LOADFISHSIZER or raises the existing
%      singleton*.
%
%      H = LOADFISHSIZER returns the handle to a new LOADFISHSIZER or the handle to
%      the existing singleton*.
%
%      LOADFISHSIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADFISHSIZER.M with the given input arguments.
%
%      LOADFISHSIZER('Property','Value',...) creates a new LOADFISHSIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loadFishSizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loadFishSizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loadFishSizer

% Last Modified by GUIDE v2.5 29-Jul-2021 10:39:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loadFishSizer_OpeningFcn, ...
                   'gui_OutputFcn',  @loadFishSizer_OutputFcn, ...
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


% --- Executes just before loadFishSizer is made visible.
function loadFishSizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loadFishSizer (see VARARGIN)

% Choose default command line output for loadFishSizer
fudgeFactor=0.55;
sz=3;
handles.zoom=0;
offset=50;
handles.tested=0;
setappdata(0,'tpath',0);
handles.cf=pwd;

rstart=getappdata(0,'restart');

if rstart==1
    fudgeFactor=getappdata(0,'ff');
    sz = getappdata(0,'dsize');
    offset = getappdata(0,'ofs')+50;
    handles.path1=getappdata(0,'pname');
    handles.calpath=getappdata(0,'cpath');
end

set(handles.thresholdLoad,'string',sprintf('%.2f',fudgeFactor));
set(handles.dilationLoad,'string',sprintf('%d',sz));
set(handles.widthLoad,'string',sprintf('%.0f',offset));

xticklabels({})
yticklabels({})

calcheck=0;
setappdata(0,'cc',calcheck);
setappdata(0,'pmm1',[]);

handles.loadbutton.String = 'Load data (uncalibrated)';
set(handles.loadbutton,'BackgroundColor','yellow');

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes loadFishSizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = loadFishSizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadtest.
function loadtest_Callback(hObject, eventdata, handles)
% hObject    handle to loadtest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fudgeFactor = str2double(get(handles.thresholdLoad,'String'));
sz = str2double(get(handles.dilationLoad,'String'));
offset = str2double(get(handles.widthLoad,'String'))-50;

oldfolder = getappdata(0,'folder1');

if offset<=-45
    offset=46;
end
if offset>=45
    offset=44;
end

rstart=getappdata(0,'restart');
if rstart==1
    oldfolder=cd(handles.path1)
end

rstart=getappdata(0,'restart');
if rstart==1
    oldfolder=cd(handles.path1)
end

[files,path1] = uigetfile( ...
{'*.jpg;*.jpeg;*.tif;*.tiff;*.bmp',...
    'Image files (*.jpg,*.jpeg,*.tif,*.tiff,*.bmp)';
   '*.*',  'All Files (*.*)'}, ...
   'Select image file for testing', ...
   'MultiSelect','off');

oldfolder = getappdata(0,'folder1');
cd(oldfolder)

if rstart==1
    cd(oldfolder)
end

if path1==0
  % user pressed cancel
  return
end
handles.tested=1;
setappdata(0,'tpath',path1);


imshow(fullfile(path1,files))
pic=imread(fullfile(path1,files));
handles.path1=path1;
handles.files=files;
handles.pic=pic;
handles.zoom=0;
guidata(hObject, handles);

% --- Executes on button press in setroibutton.
function setroibutton_Callback(hObject, eventdata, handles)
% hObject    handle to setroibutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pic=handles.pic;
cla
imshow(pic)
[xz,yz] = ginputax(handles.axes1,2);
xmin=min(xz);
xmax=max(xz);
ymin=min(yz);
ymax=max(yz);
line([xmin,xmax],[ymin,ymin],'LineWidth',2,'color','red')
line([xmin,xmax],[ymax,ymax],'LineWidth',2,'color','red')
line([xmin,xmin],[ymin,ymax],'LineWidth',2,'color','red')
line([xmax,xmax],[ymin,ymax],'LineWidth',2,'color','red')


handles.xmin=xmin;
handles.xmax=xmax;
handles.ymin=ymin;
handles.ymax=ymax;

handles.zoom=1;
guidata(hObject, handles);

% --- Executes on button press in circularroi.
function circularroi_Callback(hObject, eventdata, handles)
% hObject    handle to circularroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pic=handles.pic;
cla
imshow(pic)
[xz,yz] = ginputax(handles.axes1,2);
xmin=min(xz);
xmax=max(xz);
ymin=min(yz);
ymax=max(yz);
C=[round((xz(1)+xz(2))/2),round((yz(1)+yz(2))/2)];
X=[xz(1),yz(1);xz(2),yz(2)];
D=pdist(X,'euclidean');

handles.D=D;
handles.C=C;

viscircles(C,D/2);

handles.zoom=2;
guidata(hObject, handles);

% --- Executes on button press in deleteRoi.
function deleteRoi_Callback(hObject, eventdata, handles)
% hObject    handle to deleteRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla
pic=imread(fullfile(handles.path1,handles.files));
imshow(pic)
handles.zoom=0;
guidata(hObject, handles);

% --- Executes on button press in testparameter.
function testparameter_Callback(hObject, eventdata, handles)
% hObject    handle to testparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fudgeFactor = str2double(get(handles.thresholdLoad,'String'));
sz = str2double(get(handles.dilationLoad,'String'));
offset = str2double(get(handles.widthLoad,'String'))-50;
pic=handles.pic;
zoom=handles.zoom;
test=1;

f = waitbar(0,'Please wait...');

if zoom==1
    xmin=handles.xmin;
    xmax=handles.xmax;
    ymin=handles.ymin;
    ymax=handles.ymax;
end

if zoom==2
    C=handles.C;
    D=handles.D;
end

%% analyze data
I = rgb2gray(pic);
%I = imadjust(I,[0.2 0.8]);
Iheight=length(I(:,1));
if Iheight>=1000
    I=imgaussfilt(I,4);
    I = imsharpen(I,'radius',4,'Amount',.8);
else
       I=imgaussfilt(I,2);
    I = imsharpen(I,'radius',2,'Amount',.8); 
end

if zoom==0
[~,threshold] = edge(I(Iheight*0.2:end-Iheight*0.2,:),'sobel');
end

if zoom==1
    [~,threshold] = edge(I(round(ymin):round(ymax),round(xmin):round(xmax)),'sobel');
end

if zoom==2
    rdiff=round(cos(pi/4)*D/2);
    x1=C(2)-rdiff;
    if x1<=0
        x1=1;
    end
    x2=C(2)+rdiff
    length(I(:,1))
    if x2>=length(I(:,1))
        x2=length(I(:,1))
    end
    x3=C(1)-rdiff;
    if x3<=0
        x3=1;
    end
    x4=C(1)+rdiff;
    if x4>=length(I(1,:))
       x4=length(I(1,:))
    end
    
   [~,threshold] = edge(I(x1:x2,x3:x4),'sobel'); 
%    [~,threshold] = edge(I(C(2)-rdiff:C(2)+rdiff,C(1)-rdiff:C(1)+rdiff),'sobel');
end

BWs = edge(I,'sobel',threshold * fudgeFactor);

if zoom==1
    BWs(1:round(ymin),:)=0;
     BWs(round(ymax):end,:)=0;
     BWs(:,1:round(xmin))=0;
     BWs(:,round(xmax):end)=0;
end

if zoom==2
    [xx,yy] = meshgrid((1:size(BWs,2))-C(1),(1:size(BWs,1))-C(2)); % create x and y grid
    BWs(sqrt(xx.^2 + yy.^2) >= D/2) = 0; % set points inside the radius equal to one
end

se90 = strel('line',sz,90);
se45 = strel('line',sz,45);
se0 = strel('line',sz,0);

BWsdil = imdilate(BWs,[se90 se45 se0]);



BWdfill = imfill(BWsdil,'holes');

BWnobord = BWdfill;

seD = strel('diamond',3);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);

BWfinal =imclearborder (BWfinal,8);

Areastats = regionprops(BWfinal,'Area');
area=cell2mat(struct2cell(Areastats));
areathreshold=max(area);

if isempty(area)==0
BW2 = bwareaopen(BWfinal,round(areathreshold));

stats = regionprops('table',BW2,'Area','Orientation');
if stats.Orientation<0
    stats.Orientation=360+stats.Orientation;
end
end

close(f)

axes(handles.axes1);
if isempty('BW2')==0
imshow(labeloverlay(I,BW2,'Transparency',0.8))

if handles.zoom==1
    
    line([xmin,xmax],[ymin,ymin],'LineWidth',2,'color','red')
    line([xmin,xmax],[ymax,ymax],'LineWidth',2,'color','red')
    line([xmin,xmin],[ymin,ymax],'LineWidth',2,'color','red')
    line([xmax,xmax],[ymin,ymax],'LineWidth',2,'color','red')
end

handles.pic=pic;
end
guidata(hObject, handles);

% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fudgeFactor = str2double(get(handles.thresholdLoad,'String'));
sz = str2double(get(handles.dilationLoad,'String'));
offset = str2double(get(handles.widthLoad,'String'))-50;
zoom=handles.zoom;

if zoom==1
    pic=handles.pic;
    I = rgb2gray(pic);
    xmin=handles.xmin/length(I(1,:));
    xmax=handles.xmax/length(I(1,:));
    ymin=handles.ymin/length(I(:,1));
    ymax=handles.ymax/length(I(:,1));
end

if zoom==2
    pic=handles.pic;
    I = rgb2gray(pic);
    D=handles.D;
    C=handles.C;
end

if offset<=-45
    offset=-44;
end
if offset>=45
    offset=44;
end

setappdata(0,'ff',fudgeFactor);
setappdata(0,'dsize',sz);
setappdata(0,'ofs',offset);
setappdata(0,'z',zoom);

if zoom==1
    setappdata(0,'x1',xmin);
    setappdata(0,'x2',xmax);
    setappdata(0,'y1',ymin);
    setappdata(0,'y2',ymax);
end

if zoom==2
    setappdata(0,'diam',D);
    setappdata(0,'cent',C);
end

uiresume
guidata(hObject, handles);

close loadFishSizer

%guidata(hObject, handles);


function thresholdLoad_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdLoad as text
%        str2double(get(hObject,'String')) returns contents of thresholdLoad as a double


% --- Executes during object creation, after setting all properties.
function thresholdLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dilationLoad_Callback(hObject, eventdata, handles)
% hObject    handle to dilationLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dilationLoad as text
%        str2double(get(hObject,'String')) returns contents of dilationLoad as a double


% --- Executes during object creation, after setting all properties.
function dilationLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dilationLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function widthLoad_Callback(hObject, eventdata, handles)
% hObject    handle to widthLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthLoad as text
%        str2double(get(hObject,'String')) returns contents of widthLoad as a double


% --- Executes during object creation, after setting all properties.
function widthLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to widthLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in createcal.
function createcal_Callback(hObject, eventdata, handles)
% hObject    handle to createcal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

calibrationGUI
uiwait

calcheck = getappdata(0,'cc')

if calcheck==1
    handles.loadbutton.String = 'Load data (calibrated)';
    set(handles.loadbutton,'BackgroundColor','green');
end


guidata(hObject, handles);

% --- Executes on button press in retrievecal.
function retrievecal_Callback(hObject, eventdata, handles)
% hObject    handle to retrievecal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%currentFolder = pwd
oldfolder = getappdata(0,'folder1');
rstart=getappdata(0,'restart');
if rstart==1
    oldfolder=cd(handles.calpath)
end

testpath=getappdata(0,'tpath')

if testpath~=0
    oldfolder=cd(testpath);
end

[cfile,cpath1]=uigetfile('calibrationScale*.mat','select a calibration file');

oldfolder = getappdata(0,'folder1');
cd(oldfolder)

if cpath1==0
    cd(oldfolder)
  return
else
    
setappdata(0,'tpath',cpath1);
load(fullfile(cpath1,cfile))
setappdata(0,'cpath',cpath1);

handles.scale=scale;
if isempty(scale)==0
    setappdata(0,'pmm1',scale);
    setappdata(0,'cc',1);
end
calcheck = getappdata(0,'cc');

if calcheck==1
    handles.loadbutton.String = 'Load data (calibrated)';
    set(handles.loadbutton,'BackgroundColor','green');
end


end

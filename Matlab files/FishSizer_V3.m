function varargout = FishSizer_V3(varargin)
% FISHSIZER_V3 MATLAB code for FishSizer_V3.fig
%      FISHSIZER_V3, by itself, creates a new FISHSIZER_V3 or raises the existing
%      singleton*.
%
%      H = FISHSIZER_V3 returns the handle to a new FISHSIZER_V3 or the handle to
%      the existing singleton*.
%
%      FISHSIZER_V3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FISHSIZER_V3.M with the given input arguments.
%
%      FISHSIZER_V3('Property','Value',...) creates a new FISHSIZER_V3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FishSizer_V3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FishSizer_V3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FishSizer_V3

% Last Modified by GUIDE v2.5 26-Jul-2021 15:22:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FishSizer_V3_OpeningFcn, ...
                   'gui_OutputFcn',  @FishSizer_V3_OutputFcn, ...
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


% --- Executes just before FishSizer_V3 is made visible.
function FishSizer_V3_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FishSizer_V3 (see VARARGIN)
handles.cal=0;

folder1=pwd;

setappdata(0,'restart',0);

setappdata(0,'folder1',folder1);

set(handles.slider1,'visible','on')

Nlength=2;
Nlength2=1;

set(handles.displaylengthN,'string',sprintf('%.0f',Nlength));
set(handles.widthsegments,'string',sprintf('%.0f',Nlength2));

handles.output = hObject;

handles.folder1=folder1;
guidata(hObject, handles);

% UIWAIT makes FishSizer_V3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FishSizer_V3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete handles.datacell;

folder1=handles.folder1;

a=1;
handles.cal=0;
handles.exportbutton.String = 'Export data (uncalibrated)';
set(handles.exportbutton,'BackgroundColor','yellow');

loadFishSizer
uiwait
fudgeFactor=getappdata(0,'ff');
sz = getappdata(0,'dsize');
offset = getappdata(0,'ofs');

zoom = getappdata(0,'z');

calcheck = getappdata(0,'cc');

if calcheck==1
    handles.exportbutton.String = 'Export data (calibrated)';
    set(handles.exportbutton,'BackgroundColor','green');
    handles.scale = getappdata(0,'pmm1');
    handles.cal=1;
else
    handles.scale = 1;
end

if zoom==1
    xmin = getappdata(0,'x1');
    xmax = getappdata(0,'x2');
    ymin = getappdata(0,'y1');
    ymax = getappdata(0,'y2');
end

if zoom==2
    D = getappdata(0,'diam');
    C = getappdata(0,'cent');
end

rstart=getappdata(0,'restart');
if rstart==1
    oldfolder=cd(handles.path1);
end

testpath=getappdata(0,'tpath');
if testpath~=0
    oldfolder=cd(testpath);
end

[files,path1] = uigetfile( ...
{'*.jpg;*.jpeg;*.tif;*.tiff;*.bmp',...
    'Image files (*.jpg,*.jpeg,*.tif,*.tiff,*.bmp)';
   '*.*',  'All Files (*.*)'}, ...
   'Select image files for analysis', ...
   'MultiSelect','on');

if rstart==1 | testpath~=0
    cd(oldfolder)

end
oldfolder = getappdata(0,'folder1');
 cd(oldfolder)
 

%%
if path1~=0

tf = iscell(files);
if tf==1
    Nfiles=length(files);
    set(handles.slider1,'visible','on')
else
    Nfiles=1;
    set(handles.slider1,'visible','off')
end

if tf==1
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
end
%set(handles.slider1, 'Value', Nfiles);
set(handles.slider1, 'Value', 1);

set(handles.ThresholdD,'string',sprintf('%.2f',fudgeFactor));
%set(handles.displaythreshold,'string',sprintf('%.2f',fudgeFactor));
set(handles.szD,'string',sprintf('%d',sz));
set(handles.offsetD,'string',sprintf('%.0f',offset+50));

%%  plot
f = waitbar(0,'Please wait...');

datacell=cell(Nfiles,18);
tic
for a=Nfiles:-1:1

    v=[];
    yv=[];
if tf==1
pic=imread(fullfile(path1,cell2mat(files(a))));

filenam=cell2mat(files(a));
else
  pic=imread(fullfile(path1,files));
  filenam=files;
end


I = rgb2gray(pic);

set(handles.filenameD,'string',filenam);

Iheight=length(I(:,1));
if Iheight>=1000
    I=imgaussfilt(I,4);
    I = imsharpen(I,'radius',4,'Amount',.8);
else
       I=imgaussfilt(I,2);
    I = imsharpen(I,'radius',2,'Amount',.8); 
end

if zoom==0
    [~,threshold] = edge(I(round(Iheight*0.2):end-round(Iheight*0.2),:),'sobel');
end

if zoom==1
    I2=I(round(ymin*length(I(:,1))):round(ymax*length(I(:,1))),round(xmin*length(I(1,:))):round(xmax*length(I(1,:))));
    [~,threshold] = edge(I2,'sobel');
end

if zoom==2
    rdiff=round(cos(pi/4)*D/2);
        x1=C(2)-rdiff;
    if x1<=0
        x1=1;
    end
    x2=C(2)+rdiff;
    length(I(:,1))
    if x2>=length(I(:,1))
        x2=length(I(:,1));
    end
    x3=C(1)-rdiff;
    if x3<=0
        x3=1;
    end
    x4=C(1)+rdiff;
    if x4>=length(I(1,:))
       x4=length(I(1,:));
    end
    
   [~,threshold] = edge(I(x1:x2,x3:x4),'sobel'); 
end

BWs = edge(I,'sobel',threshold * fudgeFactor);

if zoom==1
    BWs(1:round(ymin*length(I(:,1))),:)=0;
     BWs(round(ymax*length(I(:,1))):end,:)=0;
     BWs(:,1:round(xmin*length(I(1,:))))=0;
     BWs(:,round(xmax*length(I(1,:))):end)=0;
end

if zoom==2
    [xx,yy] = meshgrid((1:size(BWs,2))-C(1),(1:size(BWs,1))-C(2)); % create x and y grid
    BWs(sqrt(xx.^2 + yy.^2) >= D/2) = 0; % set points outside the radius equal to zero
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

 stats = regionprops('table',BW2,'Area',...
     'MajorAxisLength','MinorAxisLength','Centroid','Orientation');

if stats.Orientation<0
    stats.Orientation=360+stats.Orientation;
end

Irotate=imrotate(I,-stats.Orientation);
BW2rotate=imrotate(BW2,-stats.Orientation);

[y,x]=find(BW2rotate);
xrange=linspace(min(x),max(x),max(x)-min(x)+1);
yrange=linspace(min(y),max(y),max(y)-min(y)+1);



v(:,1)=x;
v(:,2)=y;

y=length(BW2rotate(:,1))-y;


xmid=min(xrange)+round((max(x)-min(x)+1)/2);
leftidx=find(x<xmid);
leftx=x(leftidx);
lefty=y(leftidx);


rightidx=find(x>xmid);

if length(leftidx)>length(rightidx)
    or='left';
    xsep=min(xrange)+round((max(x)-min(x)+1)/4);
    tailidx=find(x>xsep);
    xtail=x(tailidx);
    ytail=y(tailidx);
    xtailrange=[xsep:max(xrange)];
    
    xhalf=x(rightidx);
    yhalf=y(rightidx);
    xhalfrange=[xmid:max(xrange)];
    phalf = polyfit(xhalf,yhalf,2);
    fhalf = polyval(phalf,xhalf,xhalfrange);
    xmid=xmid+round(length(xrange)*(offset/100));
    
else
    or='right';
    xsep=min(xrange)+round((3*(max(x)-min(x))+1)/4);
    tailidx=find(x<xsep);
    xtail=x(tailidx);
    ytail=y(tailidx);
    xtailrange=[min(xrange):xsep];

    xhalf=x(leftidx);
    yhalf=y(leftidx);
    xhalfrange=[min(xrange):xmid];
    phalf = polyfit(xhalf,yhalf,2);
    fhalf = polyval(phalf,xhalf,xhalfrange);
    xmid=xmid-round(length(xrange)*(offset/100));
    
    
end

if phalf(2)^2>0.5^2 || phalf(1)^2>0.1^2
    nfit=3;
else
    nfit=2;
end


p = polyfit(x,y,nfit);
f1=polyval(p,xrange);

CLF = hypot(diff(xrange), diff(f1));    % Calculate integrand from x,y derivatives
%CL = trapz(CLF);                   % Integrate to calculate arc length
CL = sum(CLF);

yup=zeros(1,length(xrange));
ydown=zeros(1,length(xrange));

for aa=min(xrange):max(xrange)
    if ismember(aa,x)==1
    yup(aa-min(xrange)+1)=max(y(x==aa));
    ydown(aa-min(xrange)+1)=min(y(x==aa));
    else
        yup(aa-min(xrange)+1)=NaN;
        ydown(aa-min(xrange)+1)=NaN;
    end
end

yupmedian=median(yup(round((xmid-min(xrange))-(length(xrange)/20)):round((xmid-min(xrange))+(length(xrange)/20))));
ydownmedian=median(ydown(round((xmid-min(xrange))-(length(xrange)/20)):round((xmid-min(xrange))+(length(xrange)/20))));
widthmid=yupmedian-ydownmedian;


datacell{a,1}=filenam;
datacell{a,2}=I;
datacell{a,3}=BW2;
datacell{a,4}=xrange;
datacell{a,5}=f1;
datacell{a,6}=x;
datacell{a,7}=y;
datacell{a,8}=CL;
datacell{a,9}=CL;
datacell{a,10}=xmid;
datacell{a,11}=yupmedian;
datacell{a,12}=ydownmedian;
datacell{a,13}=widthmid;
datacell{a,14}=widthmid;

%save('test4.mat')

waitbar((Nfiles-a)/Nfiles,f,'Processing your data');

else
    %msgbox('larva not detected')
    BW2=ones(length(I(:,1)),length(I(1,:)));
    widthmid=NaN;
    xmid=NaN;
    yupmedian=NaN;
    ydownmedian=NaN;
    or='left';
    
    datacell{a,1}=filenam;
    datacell{a,2}=I;
    datacell{a,3}=logical(BW2);
    datacell{a,4}=NaN;
    datacell{a,5}=NaN;
    datacell{a,6}=NaN;
    datacell{a,7}=NaN;
    datacell{a,8}=NaN;
datacell{a,9}=NaN;
datacell{a,10}=xmid;
datacell{a,11}=yupmedian;
datacell{a,12}=ydownmedian;
datacell{a,13}=widthmid;
datacell{a,14}=widthmid;
    
end





end
toc
%save('test2.mat')

close(f)

axes(handles.axes3)
cla

if exist('x','var')==1
imshow(labeloverlay(I,BW2,'Transparency',0.8))
ylim([0 length(I(:,1))])
xlim([0 length(I(1,:))])

axes(handles.axes1)
cla
axis equal
hold on
scatter(x,y,'.')
plot(xrange,f1,'LineWidth',2)
hold on
pl=line([xmid xmid],[yupmedian ydownmedian],'LineWidth',2);
pl.Color = 'green';
%width=max(y(widthmax))-min(y(widthmax));
else
    x=NaN;
    y=NaN;
    CL=NaN;
    width=NaN;
end


axes(handles.axes4)
cla

scatter(linspace(1,Nfiles,Nfiles),cell2mat(datacell(:,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
if cell2mat(datacell(a,9))==cell2mat(datacell(a,8))
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xlim([0 Nfiles+1])
ylim([0 inf])
xticks(linspace(1,length(files),Nfiles))

axes(handles.axes5)
cla

yy=datacell(a,13);
scatter(linspace(1,Nfiles,Nfiles),cell2mat(datacell(:,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xlim([0 Nfiles+1])
ylim([0 inf])
xticks(linspace(1,length(files),Nfiles))

set(handles.lengthD,'string',sprintf('%.2f',CL*handles.scale));
set(handles.widthD,'string',sprintf('%.2f',widthmid*handles.scale));

setappdata(0,'restart',1);
setappdata(0,'pname',path1);

handles.datacell=datacell;
handles.files=files;
handles.tf=tf;
handles.Nfiles=Nfiles;
handles.filenam=filenam;
handles.output = hObject;
handles.path1=path1;
handles.or=or;
handles.fudgeFactor=fudgeFactor;
handles.sz=sz;
handles.offset=offset;

end
guidata(hObject, handles);


% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacell=handles.datacell;
cal=handles.cal;
if cal==1
    scale=handles.scale;
end

% datacell{a,1}=filenam;
% datacell{a,2}=I;
% datacell{a,3}=BW2;
% datacell{a,4}=xrange;
% datacell{a,5}=f1;
% datacell{a,6}=x;
% datacell{a,7}=y;
% datacell{a,8}=CL;
% datacell{a,9}=CL;
% datacell{a,10}=xmid;
% datacell{a,11}=yupmedian;
% datacell{a,12}=ydownmedian;
% datacell{a,13}=widthmid;
% datacell{a,14}=widthmid;

T = cell2table(datacell);
T.Properties.VariableNames{1} = 'FileName';
T.Properties.VariableNames{2} = 'Image';
T.Properties.VariableNames{3} = 'BlackWhiteMask';
T.Properties.VariableNames{4} = 'Xrange';
T.Properties.VariableNames{5} = 'RegressionLine';
T.Properties.VariableNames{6} = 'XvaluesForSegmentation';
T.Properties.VariableNames{7} = 'YvaluesForSegmentation';
T.Properties.VariableNames{8} = 'Length_px';
T.Properties.VariableNames{9} = 'LengthAutomated_px';
T.Properties.VariableNames{10} = 'XvalueForWidth';
T.Properties.VariableNames{11} = 'MedianYvalueUp';
T.Properties.VariableNames{12} = 'MedianYvalueDown';
T.Properties.VariableNames{13} = 'Depth_px';
T.Properties.VariableNames{14} = 'DepthAutomated_px';

T.Length_px=round(T.Length_px);
T.LengthAutomated_px=round(T.LengthAutomated_px);
T.Depth_px=round(T.Depth_px);
T.DepthAutomated_px=round(T.DepthAutomated_px);



T1=T(:,1);
T2=T(:,8);
T3=T(:,9);
T4=T(:,13);
T5=T(:,14);

 Ttotal=[T1,T2,T3,T4,T5];

if cal==1
    
    T6=Ttotal(:,2:5);
    T6 = varfun(@(x)x.*scale, T6);
    T6.Properties.VariableNames{1} = 'length_mm';
    T6.Properties.VariableNames{2} = 'LengthAutomated_mm';
    T6.Properties.VariableNames{3} = 'Depth_mm';
    T6.Properties.VariableNames{4} = 'DepthAutomated_mm';

sd=-round(log10(scale));
        
    
    T6.length_mm=round(T6.length_mm,sd);
    T6.LengthAutomated_mm=round(T6.LengthAutomated_mm,sd);
    T6.Depth_mm=round(T6.Depth_mm,sd);
    T6.DepthAutomated_mm=round(T6.DepthAutomated_mm,sd);
    

    Ttotal=[T1,T2,T3,T4,T5,T6];
end

currentFolder = pwd;

rstart=getappdata(0,'restart');
if rstart==1
    handles.calpath = getappdata(0,'cpath');
    oldfolder=cd(handles.calpath);
end

testpath=getappdata(0,'tpath');

if testpath~=0
    oldfolder=cd(testpath);
end

newname=sprintf('DataSummary_threshold %0.2f_dill_%0.0f_offset_%0.0f.xlsx',handles.fudgeFactor,handles.sz,handles.offset);
[file,path] = uiputfile(newname);

filename = 'ManuallyChecked.xlsx';
writetable(Ttotal,fullfile(path,file))
%save('test6.mat')
if path==0
    cd(currentFolder)
  return
end

if rstart==1 | testpath~=0
    cd(oldfolder)
end
oldfolder = getappdata(0,'folder1');
cd(oldfolder)

guidata(hObject, handles);
handles.output = hObject;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
datacell=handles.datacell;
fudgeFactor=handles.fudgeFactor;
sz=handles.sz;

set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
a= round(get(handles.slider1,'Value'));

filenam=cell2mat(datacell(a,1));
I=cell2mat(datacell(a,2));
BW2=cell2mat(datacell(a,3));
xrange=cell2mat(datacell(a,4));
f1=cell2mat(datacell(a,5));
x=cell2mat(datacell(a,6));
y=cell2mat(datacell(a,7));
CL=cell2mat(datacell(a,8));
xmid=cell2mat(datacell(a,10));
yupmedian=cell2mat(datacell(a,11));
ydownmedian=cell2mat(datacell(a,12));
widthmid=cell2mat(datacell(a,13));

set(handles.filenameD,'string',filenam);


%%  plot

axes(handles.axes3)

cla

imshow(labeloverlay(I,BW2,'Transparency',0.8))
     hold on
     line(cell2mat(datacell(a,15)),cell2mat(datacell(a,16)),'LineWidth',2,'color','red')
     line(cell2mat(datacell(a,17)),cell2mat(datacell(a,18)),'LineWidth',2,'color','green')
     xlim([0,length(I(1,:))])
     ylim([0,length(I(:,1))])


axes(handles.axes1)
cla
axis equal
hold on
scatter(x,y,'.')
plot(xrange,f1,'LineWidth',2)
hold on
pl=line([xmid xmid],[yupmedian ydownmedian],'LineWidth',2);
pl.Color = 'green';

axes(handles.axes4)
cla
difflength=find(cell2mat(datacell(:,8))-cell2mat(datacell(:,9)));
eqlength=find(cell2mat(datacell(:,8))==cell2mat(datacell(:,9)));

scatter(eqlength,cell2mat(datacell(eqlength,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(difflength,cell2mat(datacell(difflength,8))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on

if cell2mat(datacell(a,9))==cell2mat(datacell(a,8))
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xticks(linspace(1,Nfiles,Nfiles))
    xlim([0 length(files)+1])
    ylim([0 inf])

else
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xticks(linspace(1,Nfiles,Nfiles))
    xlim([0 length(files)+1])
    ylim([0 inf])

end

axes(handles.axes5)
cla
diffwidth=find(cell2mat(datacell(:,13))-cell2mat(datacell(:,14)));
eqwidth=find(cell2mat(datacell(:,13))==cell2mat(datacell(:,14)));

scatter(eqwidth,cell2mat(datacell(eqwidth,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(diffwidth,cell2mat(datacell(diffwidth,13))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])

hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xlim([0 Nfiles+1])
    ylim([0 inf])
    xticks(linspace(1,length(files),Nfiles))
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xlim([0 Nfiles+1])
    ylim([0 inf])
    xticks(linspace(1,length(files),Nfiles))
    
end

set(handles.lengthD,'string',sprintf('%.2f',CL*handles.scale));
set(handles.widthD,'string',sprintf('%.2f',widthmid*handles.scale));

handles.Nfiles=Nfiles;
handles.filenam=filenam;
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in lengthbutton.
function lengthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to lengthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;

if Nfiles>1
set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
a= round(get(handles.slider1,'Value'));
else
    a=1;
end

lengthN = str2double(get(handles.displaylengthN,'String'));

[xl,yl] = ginputax(handles.axes3,lengthN+1);

datacell=handles.datacell;
I=cell2mat(datacell(a,2));
BW2=cell2mat(datacell(a,3));
axes(handles.axes3)

L = get(gca,{'xlim','ylim'});  % Get axes limits.
imshow(labeloverlay(I,BW2,'Transparency',0.8))
line(xl,yl,'LineWidth',2,'color','red')
line(cell2mat(datacell(a,17)),cell2mat(datacell(a,18)),'LineWidth',2,'color','green')
 zoom reset
set(gca,{'xlim','ylim'},L)

line(xl,yl,'LineWidth',2,'color','red')

CLF = hypot(diff(xl), diff(yl));
CL = sum(CLF);
set(handles.lengthD,'string',sprintf('%.2f',CL*handles.scale));

datacell{a,8}=CL;
datacell{a,15}=xl;
datacell{a,16}=yl;

axes(handles.axes4)
cla

difflength=find(cell2mat(datacell(:,8))-cell2mat(datacell(:,9)));
eqlength=find(cell2mat(datacell(:,8))==cell2mat(datacell(:,9)));

scatter(eqlength,cell2mat(datacell(eqlength,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(difflength,cell2mat(datacell(difflength,8))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on
scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
xticks(linspace(1,Nfiles,Nfiles))
xlim([0 length(files)+1])
ylim([0 inf])

handles.datacell=datacell;

guidata(hObject, handles);


% --- Executes on button press in clearbutton.
function clearbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
datacell=handles.datacell;


if Nfiles>1
    set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
    set(handles.slider1, 'Min', 1);
    set(handles.slider1, 'Max', Nfiles);
    a= round(get(handles.slider1,'Value'));
else
    a=1;
end

I=cell2mat(datacell(a,2));
BW2=cell2mat(datacell(a,3));
datacell{a,8}=cell2mat(datacell(a,9));
datacell{a,15}=NaN;
datacell{a,16}=NaN;

axes(handles.axes3)
L = get(gca,{'xlim','ylim'});  % Get axes limits.
imshow(labeloverlay(I,BW2,'Transparency',0.8))
line(cell2mat(datacell(a,17)),cell2mat(datacell(a,18)),'LineWidth',2,'color','green')
zoom reset
set(gca,{'xlim','ylim'},L)


axes(handles.axes4)
cla
difflength=find(cell2mat(datacell(:,8))-cell2mat(datacell(:,9)));
eqlength=find(cell2mat(datacell(:,8))==cell2mat(datacell(:,9)));

scatter(eqlength,cell2mat(datacell(eqlength,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(difflength,cell2mat(datacell(difflength,8))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on
scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
xticks(linspace(1,Nfiles,Nfiles))
xlim([0 length(files)+1])
ylim([0 inf])

set(handles.lengthD,'string',sprintf('%.2f',cell2mat(datacell(a,9))*handles.scale));


handles.datacell=datacell;


guidata(hObject, handles);


% --- Executes on button press in widthbutton.
function widthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to widthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
%datacell=handles.datacell;

%save('test.mat')

if Nfiles>1
set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
a= round(get(handles.slider1,'Value'));
else
    a=1;
end

 lengthN2 = str2num(get(handles.widthsegments,'String'));

[xl,yl] = ginputax(handles.axes3,lengthN2+1);

datacell=handles.datacell;

I=cell2mat(datacell(a,2));
BW2=cell2mat(datacell(a,3));
axes(handles.axes3)
L = get(gca,{'xlim','ylim'});  % Get axes limits.
imshow(labeloverlay(I,BW2,'Transparency',0.8))
line(cell2mat(datacell(a,15)),cell2mat(datacell(a,16)),'LineWidth',2,'color','red')
line(xl,yl,'LineWidth',2,'color','green')


zoom reset
set(gca,{'xlim','ylim'},L)

CLF2 = hypot(diff(xl), diff(yl));
CL2 = sum(CLF2);
set(handles.widthD,'string',sprintf('%.2f',CL2*handles.scale));


datacell{a,13}=CL2;
datacell{a,17}=xl;
datacell{a,18}=yl;

axes(handles.axes5)
cla
diffwidth=find(cell2mat(datacell(:,13))-cell2mat(datacell(:,14)));
eqwidth=find(cell2mat(datacell(:,13))==cell2mat(datacell(:,14)));

scatter(eqwidth,cell2mat(datacell(eqwidth,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(diffwidth,cell2mat(datacell(diffwidth,13))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])

hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xlim([0 Nfiles+1])
ylim([0 inf])
xticks(linspace(1,length(files),Nfiles))


handles.datacell=datacell;

guidata(hObject, handles);

% --- Executes on button press in clearwidthbutton.
function clearwidthbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearwidthbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
datacell=handles.datacell;

if Nfiles>1
    set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
    set(handles.slider1, 'Min', 1);
    set(handles.slider1, 'Max', Nfiles);
    a= round(get(handles.slider1,'Value'));
else
    a=1;
end

I=cell2mat(datacell(a,2));
BW2=cell2mat(datacell(a,3));

    axes(handles.axes3)
    cla
    L = get(gca,{'xlim','ylim'});  % Get axes limits.
    imshow(labeloverlay(I,BW2,'Transparency',0.8))
    line(cell2mat(datacell(a,15)),cell2mat(datacell(a,16)),'LineWidth',2,'color','red')
    zoom reset
    set(gca,{'xlim','ylim'},L)
    
datacell{a,13}=cell2mat(datacell(a,14));
datacell{a,17}=NaN;
datacell{a,18}=NaN;

axes(handles.axes5)
cla
diffwidth=find(cell2mat(datacell(:,13))-cell2mat(datacell(:,14)));
eqwidth=find(cell2mat(datacell(:,13))==cell2mat(datacell(:,14)));

hold on
scatter(diffwidth,cell2mat(datacell(diffwidth,13))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])

scatter(linspace(1,Nfiles,Nfiles),cell2mat(datacell(:,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xlim([0 Nfiles+1])
    ylim([0 inf])
    xticks(linspace(1,length(files),Nfiles))
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
    xlim([0 Nfiles+1])
    ylim([0 inf])
    xticks(linspace(1,length(files),Nfiles))
    I=cell2mat(datacell(a,2));
    BW2=cell2mat(datacell(a,3));

end

set(handles.lengthD,'string',sprintf('%.2f',cell2mat(datacell(a,9))*handles.scale));
set(handles.widthD,'string',sprintf('%.2f',cell2mat(datacell(a,13))*handles.scale));

handles.datacell=datacell;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function displayfile_Callback(hObject, eventdata, handles)
% hObject    handle to displayfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayfile as text
%        str2double(get(hObject,'String')) returns contents of displayfile as a double


% --- Executes during object creation, after setting all properties.
function displayfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displaylength_Callback(hObject, eventdata, handles)
% hObject    handle to displaylength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displaylength as text
%        str2double(get(hObject,'String')) returns contents of displaylength as a double


% --- Executes during object creation, after setting all properties.
function displaylength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaylength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displaywidth_Callback(hObject, eventdata, handles)
% hObject    handle to displaywidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displaywidth as text
%        str2double(get(hObject,'String')) returns contents of displaywidth as a double


% --- Executes during object creation, after setting all properties.
function displaywidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaywidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displaythreshold_Callback(hObject, eventdata, handles)
% hObject    handle to displaythreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displaythreshold as text
%        str2double(get(hObject,'String')) returns contents of displaythreshold as a double


% --- Executes during object creation, after setting all properties.
function displaythreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaythreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displaysize_Callback(hObject, eventdata, handles)
% hObject    handle to displaysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displaysize as text
%        str2double(get(hObject,'String')) returns contents of displaysize as a double


% --- Executes during object creation, after setting all properties.
function displaysize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displaylengthN_Callback(hObject, eventdata, handles)
% hObject    handle to displaylengthN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displaylengthN as text
%        str2double(get(hObject,'String')) returns contents of displaylengthN as a double


% --- Executes during object creation, after setting all properties.
function displaylengthN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaylengthN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function widthsegments_Callback(hObject, eventdata, handles)
% hObject    handle to widthsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthsegments as text
%        str2double(get(hObject,'String')) returns contents of widthsegments as a double


% --- Executes during object creation, after setting all properties.
function widthsegments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to widthsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function displayoffset_Callback(hObject, eventdata, handles)
% hObject    handle to displayoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayoffset as text
%        str2double(get(hObject,'String')) returns contents of displayoffset as a double


% --- Executes during object creation, after setting all properties.
function displayoffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in createbutton.
% function createbutton_Callback(hObject, eventdata, handles)
% % hObject    handle to createbutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% calibrationGUI
% guidata(hObject, handles);




% % --- Executes on button press in retrievebutton.
% function retrievebutton_Callback(hObject, eventdata, handles)
% % hObject    handle to retrievebutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% uiopen('calibrationScale*.mat')
% handles.scale=scale;
% handles.cal=1;
% handles.exportbutton.String = 'Export data (calibrated)';
% set(handles.exportbutton,'BackgroundColor','green');
% 
% guidata(hObject, handles);


% --- Executes on button press in ignorebutton.
function ignorebutton_Callback(hObject, eventdata, handles)
% hObject    handle to ignorebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
datacell=handles.datacell;
cal=handles.cal;

if Nfiles>1
set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
a= round(get(handles.slider1,'Value'));
else
    a=1;
end

if cal==1
    datacell{a,8}=NaN;
end
datacell{a,8}=NaN;
datacell{a,13}=NaN;

% save('test7.mat','datacell')

axes(handles.axes4)
cla
difflength=find(cell2mat(datacell(:,8))-cell2mat(datacell(:,9)));
eqlength=find(cell2mat(datacell(:,8))==cell2mat(datacell(:,9)));

scatter(eqlength,cell2mat(datacell(eqlength,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(difflength,cell2mat(datacell(difflength,8))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on
if cell2mat(datacell(a,9))==cell2mat(datacell(a,8))
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xticks(linspace(1,Nfiles,Nfiles))
xlim([0 length(files)+1])
ylim([0 inf])

axes(handles.axes5)
cla
diffwidth=find(cell2mat(datacell(:,13))-cell2mat(datacell(:,14)));
eqwidth=find(cell2mat(datacell(:,13))==cell2mat(datacell(:,14)));

scatter(eqwidth,cell2mat(datacell(eqwidth,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(diffwidth,cell2mat(datacell(diffwidth,13))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])

hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xlim([0 Nfiles+1])
ylim([0 inf])
xticks(linspace(1,length(files),Nfiles))

set(handles.lengthD,'string',sprintf('%.2f',cell2mat(datacell(a,8))*handles.scale));
set(handles.widthD,'string',sprintf('%.2f',cell2mat(datacell(a,13))*handles.scale));

handles.datacell=datacell;
guidata(hObject, handles);


% --- Executes on button press in includebutton.
function includebutton_Callback(hObject, eventdata, handles)
% hObject    handle to includebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nfiles=handles.Nfiles;
tf=handles.tf;
files=handles.files;
path1=handles.path1;
datacell=handles.datacell;
cal=handles.cal;

if Nfiles>1
set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', Nfiles);
a= round(get(handles.slider1,'Value'));
else
    a=1;
end

datacell{a,8}=datacell{a,9};
datacell{a,13}=datacell{a,14};

axes(handles.axes4)
cla
difflength=find(cell2mat(datacell(:,8))-cell2mat(datacell(:,9)));
eqlength=find(cell2mat(datacell(:,8))==cell2mat(datacell(:,9)));

scatter(eqlength,cell2mat(datacell(eqlength,8))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(difflength,cell2mat(datacell(difflength,8))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])
hold on
if cell2mat(datacell(a,9))==cell2mat(datacell(a,8))
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,8))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xticks(linspace(1,Nfiles,Nfiles))
xlim([0 length(files)+1])
ylim([0 inf])

axes(handles.axes5)
cla
diffwidth=find(cell2mat(datacell(:,13))-cell2mat(datacell(:,14)));
eqwidth=find(cell2mat(datacell(:,13))==cell2mat(datacell(:,14)));

scatter(eqwidth,cell2mat(datacell(eqwidth,13))*handles.scale,35,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
hold on
scatter(diffwidth,cell2mat(datacell(diffwidth,13))*handles.scale,35,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 0 1])

hold on
if cell2mat(datacell(a,14))==cell2mat(datacell(a,13))
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
else
    scatter(a,cell2mat(datacell(a,13))*handles.scale,45,'d','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[1 .2 .2])
end
xlim([0 Nfiles+1])
ylim([0 inf])
xticks(linspace(1,length(files),Nfiles))


set(handles.lengthD,'string',sprintf('%.2f',cell2mat(datacell(a,8))*handles.scale));
set(handles.widthD,'string',sprintf('%.2f',cell2mat(datacell(a,13))*handles.scale));


handles.datacell=datacell;
guidata(hObject, handles);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
  case 'a'
      Nfiles=handles.Nfiles;
      if Nfiles>1
        set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
        set(handles.slider1, 'Min', 1);
        set(handles.slider1, 'Max', Nfiles);
        a= round(get(handles.slider1,'Value'));
    else
        a=1;
      end
    if a>1
     set(handles.slider1,'Value',a-1);
    end
     slider1_Callback(hObject, eventdata, handles)
  case 's'
      Nfiles=handles.Nfiles;
      if Nfiles>1
        set(handles.slider1, 'SliderStep', [1/(Nfiles-1) , 10/(Nfiles-1) ]);
        set(handles.slider1, 'Min', 1);
        set(handles.slider1, 'Max', Nfiles);
        a= round(get(handles.slider1,'Value'));
    else
        a=1;
      end
    if a<Nfiles
      set(handles.slider1,'Value',a+1);
    end
a=a;
slider1_Callback(hObject, eventdata, handles)
  ...
% app.Slider1.Value = value; % set the slider value
% SliderValueChanging(app, event); % execute the slider callback
    
end
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

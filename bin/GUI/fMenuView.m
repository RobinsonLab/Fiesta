function fMenuView(func,varargin)
switch func
    case 'View'
        View(varargin{1},varargin{2});
    case 'ViewCheck'
        ViewCheck;
    case 'RedGreenOverlay'
        RedGreenOverlay;
    case 'DriftCorrect'
        DriftCorrect(varargin{1});
end

function View(hMainGui,idx)
if ~isempty(idx)
   hMainGui.Values.FrameIdx=idx;
else
   hMainGui.Values.FrameIdx=round(get(hMainGui.MidPanel.sFrame,'Value')); 
end
setappdata(0,'hMainGui',hMainGui);
fShow('Image');

function ViewCheck
if strcmp(get(gcbo,'Checked'),'on')==1
    set(gcbo,'Checked','off');
else
    set(gcbo,'Checked','on');
end
fShow('Image');

function RedGreenOverlay
if strcmp(get(gcbo,'Checked'),'on')==1
    set(gcbo,'Checked','off');
else
    set(gcbo,'Checked','on');
end
fShow('Image');
fShow('Tracks');

function DriftCorrect(hMainGui)
global Stack;
Drift=getappdata(hMainGui.fig,'Drift');
if ~isempty(Drift)
    if strcmp(get(hMainGui.Menu.mProjDriftCorrect,'Checked'),'off')
        d=1;
        set(hMainGui.Menu.mProjDriftCorrect,'Checked','on');
    else
        d=0;
        set(hMainGui.Menu.mProjDriftCorrect,'Checked','off');
    end
    N=length(Stack);
    MaxImage=Stack{1};
    [XG,YG]=meshgrid(1:size(MaxImage,2),1:size(MaxImage,1));
    AverageImage=double(Stack{1})*1/N;
    progressdlg('String','Drift correction of projections','Min',1,'Max',length(Stack),'Parent',hMainGui.fig);
    for n=2:length(Stack);
        Image=double(Stack{n});
        if d
            [~,m]=min(abs(Drift(:,1)-n));
            Image=interp2(XG,YG,Image,XG+Drift(m,2)/hMainGui.Values.PixSize,YG+Drift(m,3)/hMainGui.Values.PixSize);
        end
        MaxImage(:,:,2)=Image;
        MaxImage(:,:,1)=max(MaxImage,[],3);
        AverageImage=AverageImage+double(Image)*1/N;
        progressdlg(n);
    end
    MaxImage(:,:,2)=[];
    MaxImage(isnan(MaxImage))=0;
    AverageImage(isnan(AverageImage))=0;
    setappdata(hMainGui.fig,'MaxImage',MaxImage);
    setappdata(hMainGui.fig,'AverageImage',AverageImage);
    fShow('Image');
end
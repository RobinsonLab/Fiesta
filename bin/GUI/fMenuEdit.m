function fMenuEdit(func,varargin)
switch func
    case 'Find'
        Find(varargin{1});
    case 'FindNext'
        FindNext(varargin{1});     
    case 'FindMoving'
        FindMoving(varargin{1});
    case 'FindDrift'
        FindDrift(varargin{1});        
    case 'Normalize'
        Normalize(varargin{1});        
    case 'Filter'
        Filter;   
    case 'ManualTracking'
        ManualTracking(varargin{1});
    case 'ReconnectStatic'
        ReconnectStatic;
    case 'Undo'
        Undo(varargin{1});
end

function Undo(hMainGui)
global Molecule;
global Filament;
global KymoTrackMol;
global KymoTrackFil;
global BackUp;
Molecule = BackUp.Molecule;
Filament = BackUp.Filament;
KymoTrackMol = BackUp.KymoTrackMol;
KymoTrackFil = BackUp.KymoTrackFil;
fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
fRightPanel('UpdateList',hMainGui.RightPanel.pData.FilList,Filament,hMainGui.RightPanel.pData.sFilList,hMainGui.Menu.ctListFil);
set(hMainGui.Menu.mUndo,'Enable','off');
BackUp = [];
fShared('UpdateMenu',hMainGui);
fShared('ReturnFocus');
fRightPanel('UpdateKymoTracks',hMainGui);
fShow('Image');
fShow('Tracks');
if ~isempty(Molecule)||~isempty(Filament)
    set(hMainGui.MidPanel.pView,'Visible','on');
    set(hMainGui.MidPanel.pNoData,'Visible','off');
    set(hMainGui.MidPanel.tNoData,'String','No Stack or Tracks present','Visible','off');      
    drawnow expose
end

function Find(hMainGui)
global Molecule;
global Filament;
hMainGui.Search.String = fInputDlg('Find what:','');
hMainGui.Search.Mol=[];
hMainGui.Search.Fil=[];
nMol=length(Molecule);
nFil=length(Filament);
nSearchMol=1;
for i=1:nMol
    k=strfind(Molecule(i).Name,hMainGui.Search.String);
    if k>0
        hMainGui.Search.Mol(nSearchMol)=i;
        nSearchMol=nSearchMol+1;
    end
end
nSearchFil=1;
for i=1:nFil
    k=strfind(Filament(i).Name,hMainGui.Search.String);
    if k>0
        hMainGui.Search.Fil(nSearchFil)=i;
        nSearchFil=nSearchFil+1;
    end
end
hMainGui.Search.MolP=0;
hMainGui.Search.FilP=0;
if nSearchMol>1
    p=nMol-6-hMainGui.Search.Mol(1);
    if p<1
        p=1;
    end
    fMainGui('SelectObject',hMainGui,'Molecule',hMainGui.Search.Mol(1),'normal');
    getappdata(0,'hMainGui');
    set(hMainGui.RightPanel.pData.sMolList,'Value',p)
    fRightPanel('DataPanel',hMainGui);
    fRightPanel('DataMoleculesPanel',hMainGui);
    hMainGui.Search.MolP=1;
else
    if nSearchFil>1
        p=nFil-6-hMainGui.Search.Fil(1);
        if p<1
            p=1;
        end
        fMainGui('SelectObject',hMainGui,'Filament',hMainGui.Search.Fil(1),'normal');     
        getappdata(0,'hMainGui');        
        set(hMainGui.RightPanel.pData.sFilList,'Value',p)
        fRightPanel('DataPanel',hMainGui);
        fRightPanel('DataFilamentsPanel',hMainGui);
        hMainGui.Search.FilP=1;
    end
end
if nSearchMol+nSearchFil>3
    set(hMainGui.Menu.mFindNext,'Enable','on');
else    
    set(hMainGui.Menu.mFindNext,'Enable','off');
end
setappdata(0,'hMainGui',hMainGui);
fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
fRightPanel('UpdateList',hMainGui.RightPanel.pData.FilList,Filament,hMainGui.RightPanel.pData.sFilList,hMainGui.Menu.ctListFil);

function FindNext(hMainGui)
global Molecule;
global Filament;
nMol=length(Molecule);
nFil=length(Filament);
nSearchMol=length(hMainGui.Search.Mol);
nSearchFil=length(hMainGui.Search.Fil);
if nSearchMol>hMainGui.Search.MolP
    p=nMol-6-hMainGui.Search.Mol(hMainGui.Search.MolP+1);
    if p<1
        p=1;
    end
    fMainGui('SelectObject',hMainGui,'Molecule',hMainGui.Search.Mol(hMainGui.Search.MolP+1),'normal');
    getappdata(0,'hMainGui');
    set(hMainGui.RightPanel.pData.sMolList,'Value',p)
    fRightPanel('DataPanel',hMainGui);
    fRightPanel('DataMoleculesPanel',hMainGui);    
    hMainGui.Search.MolP=hMainGui.Search.MolP+1;
else
    if nSearchFil>hMainGui.Search.FilP
        p=nFil-6-hMainGui.Search.Fil(hMainGui.Search.FilP+1);
        if p<1
            p=1;
        end
        fMainGui('SelectObject',hMainGui,'Filament',hMainGui.Search.Fil(hMainGui.Search.FilP+1),'normal');     
        getappdata(0,'hMainGui');            
        set(hMainGui.RightPanel.pData.sFilList,'Value',p)
        fRightPanel('DataPanel',hMainGui);
        fRightPanel('DataFilamentsPanel',hMainGui);        
        hMainGui.Search.FilP=hMainGui.Search.FilP+1;
    end
end
if hMainGui.Search.MolP+hMainGui.Search.FilP==nSearchMol+nSearchFil
    set(hMainGui.Menu.mFindNext,'Enable','off');
end
setappdata(0,'hMainGui',hMainGui);
%fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
%fRightPanel('UpdateList',hMainGui.RightPanel.pData.FilList,Filament,hMainGui.RightPanel.pData.sFilList,hMainGui.Menu.ctListFil);

function FindMoving(hMainGui)
global Molecule;
global KymoTrackMol;
global Filament;
global KymoTrackFil;
mode=get(gcbo,'UserData');
nMol=length(Molecule);
nFil=length(Filament);
nDataMol=zeros(nMol,1);
nDisMol=zeros(nMol,1);
for i=1:nMol
    nDataMol(i)=size(Molecule(i).Results,1);
    nDisMol(i)=norm([Molecule(i).Results(nDataMol(i),3)-Molecule(i).Results(1,3) Molecule(i).Results(nDataMol(i),4)-Molecule(i).Results(1,4)]);
end
nDataFil=zeros(nFil,1);
nDisFil=zeros(nFil,1);
for i=1:nFil
    nDataFil(i)=size(Filament(i).Results,1);
    nDisFil(i)=norm([Filament(i).Results(nDataFil(i),3)-Filament(i).Results(1,3) Filament(i).Results(nDataFil(i),4)-Filament(i).Results(1,4)]);
end
if strcmp(mode,'moving') 
    answer = fInputDlg({'Enter minmum distance in nm:','Minimum number of frames'},{'100',num2str(round(max([max(nDataMol) max(nDataFil)])*0.9))});
else
    answer = fInputDlg({'Enter maxium distance in nm:','Minimum number of frames'},{'100',num2str(round(max([max(nDataMol) max(nDataFil)])*0.9))});
end
if ~isempty(answer)
    Dis = str2double(answer{1});
    mFrame = str2double(answer{2});
    for i=1:nMol
        if strcmp(mode,'moving')
            if nDisMol(i)>=Dis && nDataMol(i)>=mFrame
                Molecule = fShared('SelectOne',Molecule,KymoTrackMol,i,1);
            else
                Molecule = fShared('SelectOne',Molecule,KymoTrackMol,i,0);
            end
        else
            if nDisMol(i)<=Dis && nDataMol(i)>=mFrame
                Molecule = fShared('SelectOne',Molecule,KymoTrackMol,i,1);
            else
                Molecule = fShared('SelectOne',Molecule,KymoTrackMol,i,0);
            end
        end
    end
    for i=1:nFil
        if strcmp(mode,'moving')
            if nDisFil(i)>=Dis && nDataFil(i)>=mFrame
                Filament = fShared('SelectOne',Filament,KymoTrackFil,i,1);
            else
                Filament = fShared('SelectOne',Filament,KymoTrackFil,i,0);
            end
        else
            if nDisFil(i)<=Dis && nDataFil(i)>=mFrame
                Filament = fShared('SelectOne',Filament,KymoTrackFil,i,1);
            else
                Filament = fShared('SelectOne',Filament,KymoTrackFil,i,0);
            end
        end
    end            
    fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
    fRightPanel('UpdateList',hMainGui.RightPanel.pData.FilList,Filament,hMainGui.RightPanel.pData.sFilList,hMainGui.Menu.ctListFil);
    setappdata(0,'hMainGui',hMainGui);
    fShow('Image');
end


function FindDrift(hMainGui)
global Molecule;
global KymoTrackMol;
FileName = Molecule(1).File;
nMol = length(Molecule);
minFrame = [];
maxFrame = [];
for n = 1:nMol
    minFrame = min( [minFrame min(Molecule(n).Results(:,1))] );
    maxFrame = max( [maxFrame max(Molecule(n).Results(:,1))] );
    if ~strcmp(FileName,Molecule(n).File)
        fMsgDlg('Detected molecules of different stacks','error');
        return;
    end
end
NumDriftMol = str2double(fInputDlg('Enter number of molecules:','5'));
if ~isempty(NumDriftMol)
    Frames = (minFrame:maxFrame)';
    sFrames = length(Frames);

    %find all Molecules that have been tracked in all frames
    p=1;
    DriftMol = struct(Molecule(1));
    for n = 1:nMol
        if size(Molecule(n).Results,1) == sFrames
            if all(Molecule(n).Results(:,1) == Frames)
                DriftMol(p) = Molecule(n); 
                index{p} = n;
                drift_index{p} = p;
                p = p+1;
            end
        end
    end
    if NumDriftMol>length(DriftMol)
        fMsgDlg({'Not enough molecules for drift correction','Check whether there are enough molecules','that have been tracked in every frame'},'error');
        return;
    end
    current = [];
    nMol = length(DriftMol);
    correlation = zeros(nMol)*NaN;
    for n = 1:nMol
        for m = n+1:nMol
            X = (DriftMol(n).Results(:,3)-DriftMol(n).Results(1,3)) - (DriftMol(m).Results(:,3) - DriftMol(m).Results(1,3));
            Y = (DriftMol(n).Results(:,4)-DriftMol(n).Results(1,4)) - (DriftMol(m).Results(:,4) - DriftMol(m).Results(1,4));
            correlation(n,m) = sum( X.^2 + Y.^2 );
        end
    end
    l=1;
    while length(current) < NumDriftMol
        [~,n]=min(min(correlation,[],1));
        [~,m]=min(min(correlation,[],2));
        if n>nMol
            p = n;
        else
            p = length(DriftMol)+1;
            correlation(p,:) = NaN;
            correlation(:,p) = NaN;
        end
        index{p} = [index{n} index{m}];
        drift_index{p} = [drift_index{n}  drift_index{m}];
        current = drift_index{p};
        X = zeros(sFrames,length(current));
        Y = X;
        for k = 1:length(current)
            X(:,k) = (DriftMol(current(k)).Results(:,3)-DriftMol(current(k)).Results(1,3));
            Y(:,k) = (DriftMol(current(k)).Results(:,4)-DriftMol(current(k)).Results(1,4));
        end
        DriftMol(p).Results(:,3) = mean(X,2);
        DriftMol(p).Results(:,4) = mean(Y,2);
        for k = 1:length(DriftMol)
            X = (DriftMol(k).Results(:,3)-DriftMol(k).Results(1,3)) - (DriftMol(p).Results(:,3));
            Y = (DriftMol(k).Results(:,4)-DriftMol(k).Results(1,4)) - (DriftMol(p).Results(:,4));
            if any(ismember(drift_index{k},current))
                correlation(k,p) = NaN;
            else
                correlation(k,p) = sum( X.^2 + Y.^2 );
            end
        end
        correlation(current,p) = NaN;
        correlation(m,n) = NaN;
        correlation(p,p) = NaN;
    end
    k = find([Molecule.Selected]==1);
    for n = k
        Molecule=fShared('SelectOne',Molecule,KymoTrackMol,n,0);
    end
    for n = index{p}
        Molecule=fShared('SelectOne',Molecule,KymoTrackMol,n,1);
    end
    fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
    setappdata(0,'hMainGui',hMainGui);
end

function ReconnectStatic
global Molecule;
global Filament;
global KymoTrackMol;
global KymoTrackFil;
global Objects;
hMainGui=getappdata(0,'hMainGui');
if isempty(Objects)
    fMsgDlg({'No Objects detected!','Loading Objects could improve reconnecting'},'warning');
end  
answer = fInputDlg({'Maximum tracking error minmum in nm:','Minimum number of frames'},{'100','5'});
button = fQuestDlg({'Force tracking at static location of objects'},'FIESTA - Reconnect Static Objects',{'Yes','No'},'No');
if ~isempty(answer)
    [Molecule,KymoTrackMol,ReconnectMol]=ReconnectObj(Molecule,KymoTrackMol,str2double(answer{1}));
    [Filament,KymoTrackFil,ReconnectFil]=ReconnectObj(Filament,KymoTrackFil,str2double(answer{1}));
    if ~isempty(Objects)
        [Molecule,ReTrackMol]=AddSingleFrames(Molecule,Objects,1,str2double(answer{1}),ReconnectMol);
        [Filament,ReTrackFil]=AddSingleFrames(Filament,Objects,0,str2double(answer{1}),ReconnectFil);
    end
    if strcmp(button,'Yes')
        if isempty(Objects)
            ReTrackMol=GetReTrack(Molecule,1,ReconnectMol);
            ReTrackFil=GetReTrack(Filament,0,ReconnectFil);
        end
        [Molecule,Objects]=fManualTracking(Molecule,Objects,ReTrackMol,1,str2double(answer{1}));
        [Filament,Objects]=fManualTracking(Filament,Objects,ReTrackFil,0,str2double(answer{1}));
        [Molecule,~]=AddSingleFrames(Molecule,Objects,1,str2double(answer{1}),ReconnectMol);
        [Filament,~]=AddSingleFrames(Filament,Objects,0,str2double(answer{1}),ReconnectFil);
        
    end
    fRightPanel('UpdateList',hMainGui.RightPanel.pData.MolList,Molecule,hMainGui.RightPanel.pData.sMolList,hMainGui.Menu.ctListMol);
    fRightPanel('UpdateList',hMainGui.RightPanel.pData.FilList,Filament,hMainGui.RightPanel.pData.sFilList,hMainGui.Menu.ctListFil);
    fRightPanel('UpdateKymoTracks',hMainGui);
    fShow('Image');
    fShow('Tracks');
end

function ReTrack=GetReTrack(Object,molecule,ReconnectObj)
global Stack;
ReTrack = repmat(struct('Data',{},'Idx',[]), 1, length(Stack));
if isempty(Object)
    return;
end
mX=ReconnectObj(:,1);
mY=ReconnectObj(:,2);
for n=1:length(Stack)
    p=1;
    for m=1:length(Object)
        if~any(Object(m).Results(:,1)==n)
            if molecule
                ReTrack(n).Data{p}=[mX(m) mY(m)];
            else
                [~,t]=min(abs(Object(m).Results(:,1)-n));
                ReTrack(n).Data{p}=Object(m).Data{t}(:,1:2);
            end
            ReTrack(n).Idx(p)=m;
            p=p+1;
        end
    end
end
            
function [Object,ReTrack]=AddSingleFrames(Object,Objects,molecule,maxError,ReconnectObj)
global Config;
global Stack;
ReTrack = repmat(struct('Data',{},'Idx',[]), 1, length(Stack));
if isempty(Object)
    return;
end
mX=ReconnectObj(:,1);
mY=ReconnectObj(:,2);
r=ReconnectObj(:,3);
progressdlg('String',sprintf('Adding single frames - Frames left: %d',length(Objects)),'Max',length(Objects));
for n=1:length(Objects)
    p=1;
    if ~isempty(Objects{n})
        for m=1:length(Object)
            if ~any(Object(m).Results(:,1)==n)
                len=Objects{n}.length(1,:)==0;
                kMol=find(len==molecule);
                distance = sqrt( (Objects{n}.center_x(kMol)-mX(m)).^2 + (Objects{n}.center_y(kMol)-mY(m)).^2);
                [mDis,kDis]=min(distance);
                if mDis<3*r(m)
                    X = double([Object(m).Results(:,3); Objects{n}.center_x(kMol(kDis))]);
                    Y = double([Object(m).Results(:,4); Objects{n}.center_y(kMol(kDis))]);
                    params0 = 0.5*std(X)+0.5*std(Y);
                    R = fzero(@Likelihood2D,params0,[],X-mean(X),Y-mean(Y)); 
                    if r<maxError
                        r(m) = R;
                        mX(m) = mean(X);
                        mY(m) = mean(Y);
                        if molecule
                            Object=fShared('AddDataMol',Object,Objects,m,[],n,kMol(kDis));
                            Object(m).Results(:,5) = fDis(Object(m).Results(:,3:4));
                        else
                            Object=fShared('AddDataFil',Object,Objects,m,[],n,kMol(kDis));
                            if strcmp(Config.RefPoint,'center')==1
                                Object(m).Results(:,3:4) = Object(m).PosCenter;
                            elseif strcmp(Config.RefPoint,'start')==1
                                Object(m).Results(:,3:4) = Object(m).PosStart;
                            else
                                Object(m).Results(:,3:4) = Object(m).PosEnd;
                            end
                            Object(m).Results(:,5) = fDis(Object(m).Results(:,3:4));
                        end
                    end
                else
                    if molecule
                        ReTrack(n).Data{p}=[mX(m) mY(m)];
                    else
                        [~,t]=min(abs(Object(m).Results(:,1)-n));
                        ReTrack(n).Data{p}=Object(m).Data{t}(:,1:2);
                    end
                    ReTrack(n).Idx(p)=m;
                    p=p+1;
                end
            end
        end   
    else
        for m=1:length(Object)
            if molecule
                ReTrack(n).Data{p}=[mX(m) mY(m)];
            else
                [~,t]=min(abs(Object(m).Results(:,1)-n));
                ReTrack(n).Data{p}=Object(m).Data{t}(:,1:2);
            end
            ReTrack(n).Idx(p)=m;
            p=p+1;    
        end
    end
    progressdlg(n,sprintf('Adding single frames - Frames left: %d',length(Objects)-n));
end

function [Objects,KymoTrackObj,ReconnectObj]=ReconnectObj(Objects,KymoTrackObj,maxError)
ReconnectObj = [];
if ~isempty(Objects)
    distance=ones(length(Objects),length(Objects))*Inf;
    mX = zeros(1,length(Objects));
    mY = zeros(1,length(Objects));
    r = zeros(1,length(Objects));
    for n=1:length(Objects)
        X = double( Objects(n).Results(:,3));
        Y = double( Objects(n).Results(:,4));
        mX(n)=mean(X);
        mY(n)=mean(Y);
        distance(n,1:n-1)=sqrt( (mX(1:n-1)-mX(n) ).^2 + (mY(1:n-1)-mY(n) ).^2);
        params0 = 0.5*std(X)+0.5*std(Y);
        r(n) = fzero(@Likelihood2D,params0,[],X-mean(X),Y-mean(Y));
    end
    [~,kObj]=min(min(distance,[],1));
    [~,mObj]=min(min(distance,[],2));
    progressdlg('String',sprintf('Reconnecting - Tracks left: %d',length(Objects)),'Indeterminate','on');    
    while ~isempty(kObj) && mObj>1
        if r(kObj)>maxError
            distance(kObj,:)=[];
            distance(:,kObj)=[];
            mX(kObj)=[];
            mY(kObj)=[];
            r(kObj)=[];
            Objects(kObj) = [];
        else   
            X = double( [Objects(kObj).Results(:,3); Objects(mObj).Results(:,3)]);
            Y = double( [Objects(kObj).Results(:,4); Objects(mObj).Results(:,4)]);
            params0 = 0.5*std(X)+0.5*std(Y);
            err = fzero(@Likelihood2D,params0,[],X-mean(X),Y-mean(Y));
            if err<maxError
                Objects(kObj).Results = [Objects(kObj).Results; Objects(mObj).Results]; 
                if isfield(Objects(kObj),'Data')
                    Objects(kObj).PosCenter = [Objects(kObj).PosCenter; Objects(mObj).PosCenter]; 
                    Objects(kObj).Data = [Objects(kObj).Data; Objects(mObj).Data]; 
                end
                [m,k]=ismember([kObj mObj],[KymoTrackObj.Index]);
                if any(m)
                    if ~m(1)
                        KymoTrackObj(k(1)).Track = KymoTrackObj(k(2)).Track;
                    else
                        KymoTrackObj(k(1)).Track = sortrows([KymoTrackObj(k(1)).Track; KymoTrackObj(k(2)).Track],1);
                    end
                    idx=KymoTrackObj(k(1)).Index;
                    KymoTrackObj(k(1)).Track=sortrows(KymoTrackObj(k(1)).Track,1);
                    set(0,'CurrentFigure',hMainGui.fig);
                    set(hMainGui.fig,'CurrentAxes',hMainGui.MidPanel.aKymoGraph);
                    KymoTrackObj(k(1)).PlotHandles(1,1) = line(KymoTrackObj(k(1)).KymoTrack(:,2),KymoTrackObj(k(1)).KymoTrack(:,1),'Color',Objects(idx).Color,'Visible','off');
                    KymoTrackObj(k(1)).PlotHandles(2,1) = line(KymoTrackObj(k(1)).KymoTrack(:,2),KymoTrackObj(k(1)).KymoTrack(:,1),'Color','black','Visible','off','LineStyle','-.');        
                    KymoTrackObj(k(1)).PlotHandles(3,1) = line(KymoTrackObj(k(1)).KymoTrack(:,2),KymoTrackObj(k(1)).KymoTrack(:,1),'Color','white','Visible','off','LineStyle',':');        
                    if Objects(idx).Visible
                        set(KymoTrackObj(k(1)).PlotHandles(1,1),'Visible','on');
                    end
                    if Objects(idx).Selected
                        set(KymoTrackObj(k(1)).PlotHandles(2:3,1),'Visible','on');
                    end
                    KymoTrackObj(k(2))=[];
                end
                for n=1:length(KymoTrackObj)
                    cIndex=sum(List(2:nList)<KymoTrackObj(n).Index);
                    KymoTrackObj(n).Index=KymoTrackObj(n).Index-cIndex;
                end
                distance(mObj,:)=[];
                distance(:,mObj)=[];
                mX(kObj)=mean(X);
                mY(kObj)=mean(Y);
                r(kObj)=err;
                mX(mObj)=[];
                mY(mObj)=[];
                r(mObj)=[];
                distance(kObj,(1:kObj-1))=sqrt( (mX(1:kObj-1)-mX(kObj) ).^2 + (mY(1:kObj-1)-mY(kObj) ).^2);
                distance(kObj+1:end,kObj)=sqrt( (mX(kObj+1:end)-mX(kObj) ).^2 + (mY(kObj+1:end)-mY(kObj) ).^2);
                Objects(mObj)=[];
                Objects(kObj).Results = sortrows(Objects(kObj).Results,1); 
            else
                if distance(mObj,kObj)<r(mObj)+r(kObj)
                    distance(kObj,:)=[];
                    distance(:,kObj)=[];
                    distance(mObj,:)=[];
                    distance(:,mObj)=[];
                    mX(kObj)=[];
                    mY(kObj)=[];
                    r(kObj)=[];
                    mX(mObj)=[];
                    mY(mObj)=[];
                    r(mObj)=[];
                    Objects(mObj)=[];
                    Objects(kObj)=[];
                else
                    distance(kObj,mObj)=Inf;
                    distance(mObj,kObj)=Inf;
                end
            end

        end   
        [~,kObj]=min(min(distance,[],1));
        [~,mObj]=min(min(distance,[],2));
        progressdlg(0,sprintf('Reconnecting - Tracks left: %d',length(Objects)));
    end
    progressdlg('close');
    ReconnectObj=[mX' mY' r'];
end

function [Object,Objects]=fManualTracking(Object,Objects,ReTrackObj,molecule,maxError)
global logfile;
global Stack;
global StackInfo;
global Config;
global FiestaDir;
if isempty(ReTrackObj)
    return;
end
params.threshold = -1.5*maxError/Config.PixSize;

params.bead_model=Config.Model;
params.max_beads_per_region=1;
params.scale=Config.PixSize;
params.ridge_model = 'quadratic';

if molecule
    params.find_molecules=0;
    params.find_beads=1;
else
    params.find_molecules=1;
    params.find_beads=0;
end

params.area_threshold=Config.Threshold.Area;
params.height_threshold=0;   
params.fwhm_estimate=Config.Threshold.FWHM;
if isempty(Config.BorderMargin)
    params.border_margin = 2 * Config.Threshold.FWHM / params.scale / (2*sqrt(2*log(2)));
else
    params.border_margin = Config.BorderMargin;
end

if isempty(Config.ReduceFitBox)
    params.reduce_fit_box = 1;
else
    params.reduce_fit_box = Config.ReduceFitBox;
end

params.focus_correction = Config.FilFocus;
params.min_cod=-1;

params.binary_image_processing='none';
params.background_filter=[];

params.display = 1;

params.options = optimset( 'Display', 'off','UseParallel','never');
params.options.MaxFunEvals = 100; 
params.options.MaxIter = 10;
params.options.TolFun = 0.1;
params.options.TolX = 0.1;

if ~isempty(StackInfo)
    params.creation_time_vector = (StackInfo.CreationTime-StackInfo.CreationTime(1))/1000;
    %check wether imaging was done during change of date 
    k = params.creation_time_vector<0;
    params.creation_time_vector(k) = params.creation_time_vector(k) + 24*60*60;
end

filestr = [FiestaDir.AppData 'logfile.txt'];
logfile = fopen(filestr,'w');
progressdlg('String',sprintf('Retracking - Frames done: %d',0),'Max',length(Stack));
for n = 1:length(Stack)
    params.bw_regions = zeros(size(Stack{n}));
    for m=1:length(ReTrackObj(n).Data)
        params.bw_regions(round(ReTrackObj(n).Data{m}(:,2)/params.scale),round(ReTrackObj(n).Data{m}(:,1)/params.scale))=1;
    end
    if max(max(params.bw_regions))==1
        if length(Objects)<n 
            Objects{n}=[];
        end
        if isempty(Objects{n})
            Objects{n}=ScanImage(Stack{n},params,n);
        else
            tempObjects=ScanImage(Stack{n},params,n);
            name = fieldnames(Objects{m});
            for k = 1:length(name)
                if ~strcmp(name{k},'time')
                    tempObjects.(name{k}) = [Objects{n}.(name{k}) tempObjects.(name{k})];
                end
            end
            Objects{n}=tempObjects;
        end
    end
    progressdlg(n,sprintf('Retracking - Frames done: %d',n));
end

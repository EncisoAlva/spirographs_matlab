% Given a the list of curves created by a previous function, create a video
% file animating two circles rolling over the shape (one inside and one 
% outside) and simulating to paint its path.
%
% ---- INUPUT ------------------------------------------------------------
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%   WhCtrPos1/2  Location of the wheel center at timestamps [2x?]
%  MarkerPos1/2  Location of marker at timepoints [2x?]
%MarkerAngle1/2  Angle of marker, at each timepoint, with respect to the
%                x-axis [1x?]
%     TotalTime  Duration of animation in seconds [1]
%     AfterTime  After animation is finished, the picture will stop for
%                this amount of seconds before ending the video [1]
%       VidName  Name of video file without extension [string]
%
% ---- OUTPUT ------------------------------------------------------------
% ** No explicit output. Video is saved to path. **
%
%
function MakeVideo_Npts( nPts, WhoIsCenter, WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  CurveColor, ...
  MaxDistDelta, ...
  TotalTime, AfterTime, VidName, ExtraOpts )

nCenters = size(WhoIsCenter,2);

%%
% handle optional arguments
if ~isfield(ExtraOpts, 'LineWidth')
  ExtraOpts.LineWidth = 2;
end
if ~isfield(ExtraOpts,'Ratio')
  ExpectedRatio = 1;
else
  ExpectedRatio =  ExtraOpts.Ratio;
end
if ~isfield(ExtraOpts,'Format') 
  ExtraOpts.format = 'mp4';
end
if ~isfield(ExtraOpts,'Orientation')
  ExtraOpts.Orientation = 'in';
end

%%
% time is parametrized by the path over the Bezier curve or the path
% described by the wheel center
TimeFromCurve = ceil(nPts,1);
switch ExtraOpts.TimerefCurve
  case 'Bezier'
    for p = 1:nPts
      TimeFromCurve{p} = zeros(1,size(AllBezierPos{p},2));
      TimeFromCurve{p}(2:end) = cumsum( vecnorm( diff(AllBezierPos{p},1,2), 2, 1 ) );
    end
  case 'Wheel'
    for p = 1:nPts
      TimeFromCurve{p} = zeros(1,size(AllBezierPos{p},2));
      TimeFromCurve{p}(2:end) = cumsum( vecnorm( diff(AllWhCtrPos{p},1,2), 2, 1 ) );
    end
  case 'Average'
    for p = 1:nPts
      TimeFromCurve_pre = zeros(2,size(AllBezierPos{p},2));
      TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(AllBezierPos{p},1,2), 2, 1 ) );
      TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff( AllWhCtrPos{p},1,2), 2, 1 ) );
      TimeFromCurve{p} = mean(TimeFromCurve_pre,1);
    end
end

% duration of video
for p = 1:nPts
  TimeFromCurve{p} = TimeFromCurve{p}*((TotalTime-AfterTime)/TimeFromCurve{p}(end));
end

if isfield(ExtraOpts,'PortionFirstRound')
  for p = 1:nPts
    TimeFromCurve{p} = TimeFromCurve{p}/TimeFromCurve{p}(end);
    [FirstRoundTime,idxFirstRound] = max(TimeFromCurve{p}(TimeFromCurve{p}<=1/4));
    TimeFromCurve{p}(1:idxFirstRound) = TimeFromCurve{p}(1:idxFirstRound) * (ExtraOpts.PortionFirstRound/FirstRoundTime);
    TimeFromCurve{p}(idxFirstRound:end) = ExtraOpts.PortionFirstRound + ...
      (1-ExtraOpts.PortionFirstRound)*(TimeFromCurve{p}(idxFirstRound:end)-FirstRoundTime)/(TimeFromCurve{p}(end)-FirstRoundTime) ;
    TimeFromCurve{p} = TimeFromCurve{1}*(TotalTime-AfterTime);
  end
end

% parameters
fps     = 30;
MaxTime = ceil(TimeFromCurve{1}(end)*fps)/fps;
nTimes  = MaxTime*fps;

% help
idxx = cell(nPts,1);
for p = 1:nPts
  idxx{p} = 1:size(TimeFromCurve{p},2);
end

aang = 0:(2*pi/ ceil( 2*pi/(MaxDistDelta/WheelRadius) )):(2*pi);
circ = WheelRadius*[cos(aang); sin(aang)];
aux_angles = 0:(pi/6):(2*pi);

% original figure
close all
set(0, 'DefaultFigureColor', 'k');

f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off
%
% make sure that everything fits, and the creen ratio is ok
% if the curve is outside the path, prepare space beforehand
switch ExtraOpts.Orientation
  case 'in'
    ExtraBorder = 0;
  case 'out'
    ExtraBorder = 2*WheelRadius;
end
%
x0 = min( AllBezierPos{1}(1,:)) -ExtraBorder;
xF = max( AllBezierPos{1}(1,:)) +ExtraBorder;
y0 = min( AllBezierPos{1}(2,:)) -ExtraBorder;
yF = max( AllBezierPos{1}(2,:)) +ExtraBorder;
for p = 1:nPts
  x0 = min( x0, min(AllMarkerPos{p}(1,:)) );
  xF = max( xF, max(AllMarkerPos{p}(1,:)) );
  y0 = min( y0, min(AllMarkerPos{p}(2,:)) );
  yF = max( yF, max(AllMarkerPos{p}(2,:)) );
end
%
x_ran = xF - x0;
y_ran = yF - y0;
%
if y_ran/x_ran < ExpectedRatio
  y_ran_new = x_ran*ExpectedRatio;
  y0 = y0 - (y_ran_new - y_ran)/2;
  yF = yF + (y_ran_new - y_ran)/2;
else 
  if y_ran/x_ran > ExpectedRatio
    x_ran_new = x_ran/ExpectedRatio;
    x0 = x0 - (x_ran_new - x_ran)/2;
    xF = xF + (x_ran_new - x_ran)/2;
  end
end
%
xlim([x0 xF])
ylim([y0 yF])
if ExpectedRatio == 16/9
  set(f1,'PaperPosition',[0 0 [1080 1920]*4],'PaperUnits','points');
end
%
f2 = figure('Visible','off','Name','With circle');

% video object
switch ExtraOpts.Format
  case 'avi'
    v = VideoWriter(strcat(VidName,".avi"),'Motion JPEG AVI');
  otherwise
    v = VideoWriter(strcat(VidName,".mp4"),'MPEG-4');
end
v.Quality = 100;
open(v)

% main loop
WB = waitbar(0,strcat('Generating video (',VidName,'.mp4)...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2025)');
for i = 0:nTimes
  % 
  % update progressbar
  waitbar(i/nTimes,WB);
  if getappdata(WB,'canceling')
    disp('Ended by user.')
    close(v);
    delete(WB)
    break
  end
  %
  CurrPts = cell(nPts,1);
  for p = 1:nPts
    CurrPts{p} = idxx{p}(TimeFromCurve{p}<=i/fps);
  end
  if ~isempty(CurrPts{1}) % if no points will be added. skip drawing loop
  %
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  plot(DecorativeBez(1,:),DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)
  %
  % add all strokes of the marker up to the current time
  for p = 1:nPts
    plot(AllMarkerPos{p}(1,CurrPts{p}),AllMarkerPos{p}(2,CurrPts{p}),'Color',CurveColor{p},'LineWidth',ExtraOpts.LineWidth)
  end
  %
  jj = zeros(nPts,1);
  for p = 1:nPts
    jj(p) = min( [max(CurrPts{p}), size(AllLocTime{p},2), size(AllBezierPos{p},2)]);
  end
  %
  % locate center of wheels
  RefWheelCtr = zeros(2,nCenters);
  RefAngle = zeros(1,nCenters);
  for q = 1:nCenters
    RefWheelCtr(:,q) = AllWhCtrPos{WhoIsCenter(q)}(:,jj(WhoIsCenter(q)));
    RefAngle(q) =   AllMarkerAngle{WhoIsCenter(q)}(  jj(WhoIsCenter(q)));
  end
  %
  % the white dot is unnecessarily challenging
  compTime = -Inf;
  for p = 1:nPts
    if AllLocTime{p}(jj(p)) > compTime
      compTime = AllLocTime{p}(  jj(p));
      RefBez = AllBezierPos{p}(:,jj(p));
    end
  end
  %
  for q = 1:nCenters
    fill(RefWheelCtr(1,q)+circ(1,:),RefWheelCtr(2,q)+circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
    for w = 1:size(aux_angles,2)
      plot(RefWheelCtr(1,q)+[0,cos(aux_angles(w)+RefAngle(q))*WheelRadius],...
           RefWheelCtr(2,q)+[0,sin(aux_angles(w)+RefAngle(q))*WheelRadius],...
        'Color',[0,0,0,0.5])
    end
  end
  for p = 1:nPts
    scatter(AllMarkerPos{p}(1,jj(p)),AllMarkerPos{p}(2,jj(p)),10,CurveColor{p},'filled')
  end
  scatter(RefBez(1),RefBez(2),10,'white','filled')
  %
  end
  writeVideo(v,getframe)
end

%%
% stop for some time after everything is finished, WITH the wheel on
for stopper = 0:(fps*(AfterTime/2))
  writeVideo(v,getframe)
end

% plot again, without the wheel
for tmp = 1:1
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  % add all strokes of the marker up to the current time
  for p = 1:nPts
    plot(AllMarkerPos{p}(1,CurrPts{p}),AllMarkerPos{p}(2,CurrPts{p}),'Color',CurveColor{p},'LineWidth',ExtraOpts.LineWidth)
  end
  %
end

%%
% stop for some time after everything is finished
for stopper = 0:(fps*(AfterTime/2))
  writeVideo(v,getframe)
end

% save as picture, it may be useful
[~,name_without_extension, ~] = fileparts(Vidname);
saveas(f2, strcat(name_without_extension,'.png'));

% delete waitbar, if it is still around
if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% finalize the video object and close the figure
close(v);

end
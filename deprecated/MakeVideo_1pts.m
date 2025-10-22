% Given a the list of curves created by a previous function, create a video
% file animating a circle rolling over the shape and simulating to paint
% its path.
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
function MakeVideo_1pts( WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  CurveColor, ...
  MaxDistDelta, ...
  TotalTime, AfterTime, VidName, ExtraOpts )

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
% time, parametrized by the arc of the marker or the ar of the wheel center
switch ExtraOpts.TimerefCurve
  case 'Bezier'
    TimeFromCurve = zeros(1,size(AllBezierPos{1},2));
    TimeFromCurve(2:end) = cumsum( vecnorm( diff(AllBezierPos{1},1,2), 2, 1 ) );
  case 'Wheel'
    TimeFromCurve = zeros(1,size(AllWhCtrPos{1},2));
    TimeFromCurve(2:end) = cumsum( vecnorm( diff(AllWhCtrPos{1},1,2), 2, 1 ) );
  case 'Average'
    TimeFromCurve_pre = zeros(2,size(AllBezierPos{1},2));
    TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(AllBezierPos{1},1,2), 2, 1 ) );
    TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(AllWhCtrPos{1}, 1,2), 2, 1 ) );
    TimeFromCurve = mean(TimeFromCurve_pre,1);
end

% duration of video
TimeFromCurve = TimeFromCurve*((TotalTime-AfterTime)/TimeFromCurve(end));

% parameters
fps = 30;
%MaxTime = ceil(TimeFromMarker(end)*fps)/fps;
MaxTime = ceil(TimeFromCurve(end)*fps)/fps;
nTimes = MaxTime*fps;

% help
idxx = 1:size(TimeFromCurve,2);
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
x0 = min( [min(AllMarkerPos{1}(1,:)), min(AllMarkerPos{2}(1,:)), (min(AllBezierPos{1}(1,:))-ExtraBorder)] );
xF = max( [max(AllMarkerPos{1}(1,:)), max(AllMarkerPos{2}(1,:)), (max(AllBezierPos{1}(1,:))+ExtraBorder)] );
y0 = min( [min(AllMarkerPos{1}(2,:)), min(AllMarkerPos{2}(2,:)), (min(AllBezierPos{1}(2,:))-ExtraBorder)] );
yF = max( [max(AllMarkerPos{1}(2,:)), max(AllMarkerPos{2}(2,:)), (max(AllBezierPos{1}(2,:))+ExtraBorder)] );
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

% prepare finished figure for sneek peek
% take specifications from figure 1
%copyobj(f1.Children,f2)
%set(0,"CurrentFigure",f2)
%plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),'Color',CurveColor{2})
%plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),'Color',CurveColor{1})

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
  %
  % take specifications from figure 1
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  plot(DecorativeBez(1,:),DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)
  %
  CurrPts = idxx(TimeFromCurve<=(i+0.1)/fps);
  if ~isempty(CurrPts) % if no points will be added. skip drawing loop
  %
  % add a few strokes of the marker, then copy to figure 2
  plot(AllMarkerPos{1}(1,CurrPts),AllMarkerPos{1}(2,CurrPts),'Color',CurveColor{2},'LineWidth',ExtraOpts.LineWidth)
  %
  j = max(CurrPts);
    RefWheelCtr = AllWhCtrPos{1}(:,j);
    RefAngle = AllMarkerAngle{1}(j);
    RefBez = AllBezierPos{1}(:,j);
  %
  fill(RefWheelCtr(1)+circ(1,:),RefWheelCtr(2)+circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.25); 
  for w = 1:size(aux_angles,2)
    plot(RefWheelCtr(1)+[0,cos(aux_angles(w)+RefAngle)*WheelRadius],RefWheelCtr(2)+[0,sin(aux_angles(w)+RefAngle)*WheelRadius],...
      'Color',[0,0,0,0.5])
  end
  scatter(AllMarkerPos{1}(1,j),AllMarkerPos{1}(2,j),10,CurveColor{2},'filled')
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
for i = 1:1
  % take specifications from figure 1
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  % add a few strokes of the marker, then copy to figure 2
  plot(AllMarkerPos{1}(1,CurrPts),AllMarkerPos{1}(2,CurrPts),CurveColor{2},'LineWidth',ExtraOpts.LineWidth)
end


%%
% stop for some time after everything is finished
for stopper = 0:(fps*(AfterTime/2))
  writeVideo(v,getframe)
end

if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% finalize the video object and close the figure
close(v);

end
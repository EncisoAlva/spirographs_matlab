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
function MakeVideo( obj, VidName, VideoOpts, varargin )

% determine duration of the video
if size(varargin,2) > 0
  TotalTime = varargin{1};
else
  TotalTime = 30;
end
if size(varargin,2) > 0
  AfterTime = varargin{2};
else
  AfterTime = 7.5;
end

% patch
nMarkers = 1;

% parameters
nCenters = size(VideoOpts.WhoIsCenter,2);
if size(VideoOpts.WheelRadii,2) < max(VideoOpts.WhoIsCenter)
  WhoIsCenter_mod = mod(VideoOpts.WhoIsCenter-1, size(VideoOpts.WheelRadii,2)) + 1;
  CenterRadii = VideoOpts.WheelRadii(WhoIsCenter_mod);
else
  CenterRadii = VideoOpts.WheelRadii(VideoOpts.WhoIsCenter);
end

%%
% handle optional arguments
if ~isfield(VideoOpts, 'LineWidth')
  VideoOpts.LineWidth = 2;
end
if ~isfield(VideoOpts,'Ratio')
  ExpectedRatio = 1;
else
  ExpectedRatio =  VideoOpts.Ratio;
end
if ~isfield(VideoOpts,'Format') 
  VideoOpts.format = 'mp4';
end
if ~isfield(VideoOpts,'Orientation')
  VideoOpts.Orientation = 'in';
end
if ~isfield(VideoOpts, 'BackgroundColor')
  VideoOpts.BackgroundColor = 'black';
end
if ~isfield(VideoOpts, 'FillBezier')
  VideoOpts.FillBezier      = false;
end
if ~isfield(VideoOpts, 'FillBezierColor')
  VideoOpts.FillBezierColor = 'black'; % non-empty color
end
if ~isfield(VideoOpts, 'FillMarkerCurve')
  VideoOpts.FillMarkerCurve = false;
end

%%
% time is parametrized by the path over the Bezier curve or the path
% described by the wheel center
nPts = size(obj.BezierPos,2);
switch VideoOpts.TimeRefCurve
  case 'Bezier'
    TimeFromCurve = zeros(1,nPts);
    TimeFromCurve(2:end) = cumsum( vecnorm( diff(obj.BezierPos,1,2), 2, 1 ) );
  case 'Marker'
    TimeFromCurve = zeros(1,nPts);
    TimeFromCurve(2:end) = cumsum( vecnorm( diff(obj.MarkerPos,1,2), 2, 1 ) );
  case 'Wheel'
    TimeFromCurve = zeros(1,nPts);
    TimeFromCurve(2:end) = cumsum( vecnorm( diff(obj.WhCtrPos,1,2), 2, 1 ) );
  case 'Average'
    TimeFromCurve_pre = zeros(2,nPts);
    TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(obj.BezierPos,1,2), 2, 1 ) );
    TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff( obj.WhCtrPos,1,2), 2, 1 ) );
    TimeFromCurve = mean(TimeFromCurve_pre,1);
  case 'Avg_MarkerBezier'
    TimeFromCurve_pre = zeros(2,nPts);
    TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(obj.MarkerPos,1,2), 2, 1 ) );
    TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff( obj.WhCtrPos,1,2), 2, 1 ) );
    TimeFromCurve = mean(TimeFromCurve_pre,1);
end

% duration of video
TimeFromCurve = TimeFromCurve*((TotalTime-AfterTime)/TimeFromCurve(end));

% parameters
fps     = 30;
MaxTime = ceil(TimeFromCurve(end)*fps)/fps;
nTimes  = MaxTime*fps;

% help
idxx = 1:nPts;

% rolling wheel, as a polygon
aang = 0:(2*pi/ ceil( 2*pi/(obj.Tol/max(CenterRadii)) )):(2*pi);
aux_angles = 0:(pi/6):(2*pi);
switch obj.Method
  case 'Default'
    circ = [cos(aang); sin(aang)];
  case 'Hole'
    circ = [ [cos(aang); sin(aang)], obj.DecorativeHole ];
end

% original figure
close all
set(0, 'DefaultFigureColor', VideoOpts.BackgroundColor);

f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off
%
% make sure that everything fits, and the creen ratio is ok
% if the curve is outside the path, prepare space beforehand
switch VideoOpts.Orientation
  case 'in'
    ExtraBorder = 0;
  case 'out'
    ExtraBorder = 2*max(CenterRadii);
end
%
x0 = min(obj.DecorativeBez(1,:));
xF = max(obj.DecorativeBez(1,:));
y0 = min(obj.DecorativeBez(2,:));
yF = max(obj.DecorativeBez(2,:));
%
x0 = min( x0, min(obj.MarkerPos(1,:)) );
xF = max( xF, max(obj.MarkerPos(1,:)) );
y0 = min( y0, min(obj.MarkerPos(2,:)) );
yF = max( yF, max(obj.MarkerPos(2,:)) );
%
x0 = x0 -ExtraBorder;
xF = xF +ExtraBorder;
y0 = y0 -ExtraBorder;
yF = yF +ExtraBorder;
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
    x_ran_new = y_ran/ExpectedRatio;
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
switch VideoOpts.Format
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
  CurrPts = idxx(TimeFromCurve<=i/fps);
  if ~isempty(CurrPts) % if no points will be added. skip drawing loop
  %
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  if VideoOpts.FillBezier
    fill(obj.DecorativeBez(1,:),obj.DecorativeBez(2,:), VideoOpts.FillBezierColor, 'EdgeColor', 'none')
  else 
    plot(obj.DecorativeBez(1,:),obj.DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)
  end
  %
  % add all strokes of the marker up to the current time
  if obj.Multicolor
    % add NaN so that Matlab understands that it is not a closed curve
    drawPts = [ obj.MarkerPos(:,CurrPts), NaN(2,1)];
    drawCol = [ obj.ColorFunc(:,CurrPts)'; NaN(1,3)];
    fill(drawPts(1,:),drawPts(2,:), 'k',...
      'FaceVertexCData',drawCol,'EdgeColor','interp','LineWidth',VideoOpts.LineWidth)
  else
    plot(obj.MarkerPos(1,CurrPts),obj.MarkerPos(2,CurrPts),...
      'Color',obj.ColorVector,'LineWidth',VideoOpts.LineWidth)
  end
  %
  jj = min( [max(CurrPts), nPts]);
  %
  % locate center of wheels
  RefWheelCtr = obj.WhCtrPos(:,jj);
  RefAngle    = obj.AllMarkerAngle(jj);
  %
  % the white dot is unnecessarily challenging
  RefBez = obj.BezierPos(:,jj);
  %
  for q = 1:nCenters
    switch VideoOpts.Method
      case 'Default'
        fill(RefWheelCtr(1,q)+circ(1,:)*VideoOpts.WheelRadii(q),RefWheelCtr(2,q)+circ(2,:)*VideoOpts.WheelRadii(q), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
      case 'Hole'
        th = RefAngle(q);
        rolled_circ = [cos(th), -sin(th); sin(th), cos(th)] * circ * VideoOpts.WheelRadii(q);
        fill(RefWheelCtr(1,q)+rolled_circ(1,:),RefWheelCtr(2,q)+rolled_circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
        %
        rolled_hole = [cos(th), -sin(th); sin(th), cos(th)] * VideoOpts.DecorativeHole*VideoOpts.WheelRadii(q);
        plot(RefWheelCtr(1,q)+rolled_hole(1,:), RefWheelCtr(2,q)+rolled_hole(2,:), 'Color',[0,0,0,0.5])
    end
    for w = 1:size(aux_angles,2)
      plot(RefWheelCtr(1,q)+[0,cos(aux_angles(w)+RefAngle(q))*CenterRadii(q)],...
           RefWheelCtr(2,q)+[0,sin(aux_angles(w)+RefAngle(q))*CenterRadii(q)],...
        'Color',[0,0,0,0.5])
    end
  end
  if obj.Multicolor
    scatter(obj.MarkerPos(1,jj),obj.MarkerPos(2,jj),20,obj.ColorFunc(:,jj)','filled')
  else 
    scatter(obj.MarkerPos(1,jj),obj.MarkerPos(2,jj),20,obj.ColorVector,'filled')
  end
  scatter(RefBez(1),RefBez(2),10,'white','filled')
  %
  end
  writeVideo(v,getframe)
end

%%
% stop for some time after everything is finished, WITH the wheel on
for stopper = 0:(fps*(AfterTime/4))
  writeVideo(v,getframe)
end

% plot again, without the wheel
for tmp = 1:1
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  % add bezier path ONLY if requested
  if VideoOpts.FillBezier
    fill(obj.DecorativeBez(1,:),obj.DecorativeBez(2,:), VideoOpts.FillBezierColor, 'EdgeColor', 'none')
  end
  % add all strokes of the marker up to the current time
  if VideoOpts.FillMarkerCurve
    fill(obj.MarkerPos(1,:),obj.MarkerPos(2,:),obj.ColorVector, 'EdgeColor', 'none')
  else
    if obj.Multicolor
      % add NaN so that Matlab understands that it is not a closed curve
      drawPts = [ obj.MarkerPos, NaN(2,1)];
      drawCol = [ obj.ColorFunc'; NaN(1,3)];
      fill(drawPts(1,:),drawPts(2,:), 'k',...
        'FaceVertexCData',drawCol,'EdgeColor','interp','LineWidth',VideoOpts.LineWidth)
    else
      plot(obj.MarkerPos(1,:),obj.MarkerPos(2,:),...
        'Color',obj.ColorVector,'LineWidth',VideoOpts.LineWidth)
    end
  end
  %
end

%%
% stop for some time after everything is finished
for stopper = 0:(fps*(3*AfterTime/4))
  writeVideo(v,getframe)
end

% save as picture, it may be useful
[~,name_without_extension, ~] = fileparts(VidName);
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
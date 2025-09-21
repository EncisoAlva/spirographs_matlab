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
function MakeVideo_2pts( WheelRadius, ...
  BezierPos, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2,...
  MaxDistDelta, ...
  TotalTime, AfterTime, VidName )

% time, parametrized by the arc of the marker or the ar of the wheel center
%TimeFromMarker = zeros(1,size(MarkerPos1,2));
%TimeFromMarker(2:end) = cumsum( vecnorm( diff(MarkerPos1,1,2), 2, 1 ) );

TimeFromWheel1 = zeros(1,size(WhCtrPos1,2));
TimeFromWheel1(2:end) = cumsum( vecnorm( diff(WhCtrPos1,1,2), 2, 1 ) );

TimeFromWheel2 = zeros(1,size(WhCtrPos2,2));
TimeFromWheel2(2:end) = cumsum( vecnorm( diff(WhCtrPos2,1,2), 2, 1 ) );

% duration of video
TimeFromWheel1 = TimeFromWheel1*((TotalTime-AfterTime)/TimeFromWheel1(end));
TimeFromWheel2 = TimeFromWheel2*((TotalTime-AfterTime)/TimeFromWheel2(end));

% parameters
fps = 30;
%MaxTime = ceil(TimeFromMarker(end)*fps)/fps;
MaxTime = ceil(TimeFromWheel1(end)*fps)/fps;
nTimes = MaxTime*fps;

% help
%idxx = 1:size(TimeFromMarker,2);
idxx1 = 1:size(TimeFromWheel1,2);
idxx2 = 1:size(TimeFromWheel2,2);
aang = 0:(2*pi/ ceil( 2*pi/(MaxDistDelta/WheelRadius) )):(2*pi);
circ = WheelRadius*[cos(aang); sin(aang)];
aux_angles = 0:(pi/3):(2*pi);

% original figure
close all
set(0, 'DefaultFigureColor', 'k');

f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off
xlim([min(MarkerPos1(1,:)) max(MarkerPos1(1,:))])
ylim([min(MarkerPos1(2,:)) max(MarkerPos1(2,:))])
xlim([min(MarkerPos2(1,:)) max(MarkerPos2(1,:))])
ylim([min(MarkerPos2(2,:)) max(MarkerPos2(2,:))])
%
fill(BezierPos(1,:),BezierPos(2,:),'Color',[.3 .3 .3], 'EdgeColor', 'none'); 
%
f2 = figure('Visible','off','Name','With circle');

% video object
v = VideoWriter(strcat(VidName,".mp4"),'MPEG-4');
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
  %
  % take specifications from figure 1
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  CurrPts1 = idxx1(TimeFromWheel1<=(i+0.1)/fps);
  CurrPts2 = idxx2(TimeFromWheel2<=(i+0.1)/fps);
  if ~( isempty(CurrPts1) & isempty(CurrPts2) ) % if no points will be added. skip drawing loop
  %
  % add a few strokes of the marker, then copy to figure 2
  plot(MarkerPos1(1,CurrPts1),MarkerPos1(2,CurrPts1),'magenta')
  plot(MarkerPos2(1,CurrPts2),MarkerPos2(2,CurrPts2),'yellow')
  %
  j1 = max(CurrPts1);
  j2 = max(CurrPts2);
  if ~isempty(j1)
    RefWheelCtr = WhCtrPos1(:,j1);
    RefAngle = MarkerAngle1(j1);
  else
    RefWheelCtr = WhCtrPos2(:,j2);
    RefAngle = MarkerAngle2(j2);
  end
  fill(RefWheelCtr(1)+circ(1,:),RefWheelCtr(2)+circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
  for w = 1:size(aux_angles,2)
    plot(RefWheelCtr(1)+[0,cos(aux_angles(w)+RefAngle)*WheelRadius],RefWheelCtr(2)+[0,sin(aux_angles(w)+RefAngle)*WheelRadius],...
      'Color',[0,0,0,0.5])
  end
  scatter(MarkerPos1(1,j1),MarkerPos1(2,j1),10,'magenta','filled')
  scatter(MarkerPos2(1,j2),MarkerPos2(2,j2),10,'yellow','filled')
  %
  end
  writeVideo(v,getframe)
end

if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% stop for some time after everything is finished
for stopper = 0:(fps*AfterTime)
  writeVideo(v,getframe)
end

% finalize the video object and close the figure
close(v);

end
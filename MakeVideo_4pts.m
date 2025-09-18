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
function MakeVideo_4pts( WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  MaxDistDelta, ...
  TotalTime, AfterTime, VidName )

% time is parametrized by the path over the Bezier curve
TimeFromBez1A = zeros(1,size(AllBezierPos{1},2));
TimeFromBez1A(2:end) = cumsum( vecnorm( diff(AllBezierPos{1},1,2), 2, 1 ) );

TimeFromBez2A = zeros(1,size(AllBezierPos{2},2));
TimeFromBez2A(2:end) = cumsum( vecnorm( diff(AllBezierPos{2},1,2), 2, 1 ) );

TimeFromBez1B = zeros(1,size(AllBezierPos{3},2));
TimeFromBez1B(2:end) = cumsum( vecnorm( diff(AllBezierPos{3},1,2), 2, 1 ) );

TimeFromBez2B = zeros(1,size(AllBezierPos{4},2));
TimeFromBez2B(2:end) = cumsum( vecnorm( diff(AllBezierPos{4},1,2), 2, 1 ) );

% duration of video
TimeFromBez1A = TimeFromBez1A*((TotalTime-AfterTime)/TimeFromBez1A(end));
TimeFromBez2A = TimeFromBez2A*((TotalTime-AfterTime)/TimeFromBez2A(end));

TimeFromBez1B = TimeFromBez1B*((TotalTime-AfterTime)/TimeFromBez1B(end));
TimeFromBez2B = TimeFromBez2B*((TotalTime-AfterTime)/TimeFromBez2B(end));

% parameters
fps = 30;
%MaxTime = ceil(TimeFromMarker(end)*fps)/fps;
MaxTime = ceil(TimeFromBez1A(end)*fps)/fps;
nTimes = MaxTime*fps;

% help
%idxx = 1:size(TimeFromMarker,2);
idxx1A = 1:size(TimeFromBez1A,2);
idxx2A = 1:size(TimeFromBez2A,2);
idxx1B = 1:size(TimeFromBez1B,2);
idxx2B = 1:size(TimeFromBez2B,2);

aang = 0:(2*pi/ ceil( 2*pi/(MaxDistDelta/WheelRadius) )):(2*pi);
circ = WheelRadius*[cos(aang); sin(aang)];
aux_angles = 0:(pi/3):(2*pi);

% original figure
close all
f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off
xlim([min(AllMarkerPos{1}(1,:)) max(AllMarkerPos{1}(1,:))])
ylim([min(AllMarkerPos{1}(2,:)) max(AllMarkerPos{1}(2,:))])
xlim([min(AllMarkerPos{2}(1,:)) max(AllMarkerPos{2}(1,:))])
ylim([min(AllMarkerPos{2}(2,:)) max(AllMarkerPos{2}(2,:))])
%
%fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
plot(DecorativeBez(1,:),DecorativeBez(2,:),'k','LineWidth',2)
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
  %CurrPts = idxx(and(TimeFromMarker>=(i-1.1)/fps,TimeFromMarker<=(i+0.1)/fps));
  CurrPts1A = idxx1A(TimeFromBez1A<=(i+0.1)/fps);
  CurrPts2A = idxx2A(TimeFromBez2A<=(i+0.1)/fps);
  CurrPts1B = idxx1B(TimeFromBez1B<=(i+0.1)/fps);
  CurrPts2B = idxx2B(TimeFromBez2B<=(i+0.1)/fps);
  if ~( isempty(CurrPts1A) & isempty(CurrPts2A) ) % if no points will be added. skip drawing loop
  %
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)
  %
  % add all strokes of the marker up to the current time
  plot(AllMarkerPos{1}(1,CurrPts1A),AllMarkerPos{1}(2,CurrPts1A),'magenta')
  plot(AllMarkerPos{2}(1,CurrPts2A),AllMarkerPos{2}(2,CurrPts2A),'yellow')
  plot(AllMarkerPos{3}(1,CurrPts1B),AllMarkerPos{3}(2,CurrPts1B),'yellow')
  plot(AllMarkerPos{4}(1,CurrPts2B),AllMarkerPos{4}(2,CurrPts2B),'magenta')
  %
  j1A = max(CurrPts1A);
  j2A = max(CurrPts2A);
  j1B = max(CurrPts1B);
  j2B = max(CurrPts2B);
  if ~isempty(j1A)
    RefWheelCtrA = AllWhCtrPos{1}(:,j1A);
    RefAngleA = AllMarkerAngle{1}(j1A);
  else
    RefWheelCtrA = AllWhCtrPos{2}(:,j2A);
    RefAngleA = AllMarkerAngle{2}(j2A);
  end
  if ~isempty(j1B)
    RefWheelCtrB = AllWhCtrPos{3}(:,j1B);
    RefAngleB = AllMarkerAngle{3}(j1B);
  else
    RefWheelCtrB = AllWhCtrPos{4}(:,j2B);
    RefAngleB = AllMarkerAngle{4}(j2B);
  end
  %
  % the white dot is unnecessarily challenging
  compTimes  = [LocTime1A(j1A), LocTime1B(j1B), LocTime2A(j2A), LocTime2B(j2B)];
  [~,refIdx] = max(compTimes);
  compRefBez = [AllBezierPos{1}(:,j1A), AllBezierPos{3}(:,j1B), AllBezierPos{2}(:,j2A), AllBezierPos{4}(:,j2B)];
  RefBez = compRefBez(:,refIdx);
  %
  fill(RefWheelCtrA(1)+circ(1,:),RefWheelCtrA(2)+circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
  fill(RefWheelCtrB(1)+circ(1,:),RefWheelCtrB(2)+circ(2,:), 'cyan', 'EdgeColor', 'none','FaceAlpha',0.15); 
  for w = 1:size(aux_angles,2)
    plot(RefWheelCtrA(1)+[0,cos(aux_angles(w)+RefAngleA)*WheelRadius],RefWheelCtrA(2)+[0,sin(aux_angles(w)+RefAngleA)*WheelRadius],...
      'Color',[0,0,0,0.5])
    plot(RefWheelCtrB(1)+[0,cos(aux_angles(w)+RefAngleB)*WheelRadius],RefWheelCtrB(2)+[0,sin(aux_angles(w)+RefAngleB)*WheelRadius],...
      'Color',[0,0,0,0.5])
  end
  scatter(AllMarkerPos{1}(1,j1A),AllMarkerPos{1}(2,j1A),10,'magenta','filled')
  scatter(AllMarkerPos{2}(1,j2A),AllMarkerPos{2}(2,j2A),10,'yellow','filled')
  scatter(AllMarkerPos{3}(1,j1B),AllMarkerPos{3}(2,j1B),10,'yellow','filled')
  scatter(AllMarkerPos{4}(1,j2B),AllMarkerPos{4}(2,j2B),10,'magenta','filled')
  scatter(RefBez(1),RefBez(2),10,'white','filled')
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
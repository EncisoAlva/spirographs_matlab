% the code for the curve is known to work, this should be the first real
% example
% the sketch of the moon took me some time, as I did it analytically

% designer stuff
MarkerAngle0 = pi;

% technical stuff
MaxDistDelta = 0.001;
CloseTol = 0.001;
MaxSpins = 50;

% specially crafted control points
CtrlPtsArray = {[...
  [-1,0]', ...
  [-1,0]'+(4/3)*(tan(pi/8)/sin(pi/2))*[0,1]',...
  [0,1]'-(4/3)*(tan(pi/8)/sin(pi/2))*[sin(pi/2),cos(pi/2)]',...
  [0,1]'...
  ],[...
  [0,1]', ...
  [0, 1]'-(4/3)*(tan(pi/6)/sin(pi/3))*[sin(pi/3), cos(pi/3)]',...
  [0,-1]'-(4/3)*(tan(pi/6)/sin(pi/3))*[sin(pi/3),-cos(pi/3)]',...
  [0,-1]'...
  ],[...
  [0,-1]', ...
  [0,-1]'-(4/3)*(tan(pi/8)/sin(pi/2))*[sin(pi/2),-cos(pi/2)]',...
  [-1,0]'-(4/3)*(tan(pi/8)/sin(pi/2))*[0,1]',...
  [-1,0]' ...
  ]};

% specific to this example
WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/15;
MarkerRadius = WheelRadius*(1);

% Bezier curve
BezierPos = AllBezierEval( CtrlPtsArray, MaxDistDelta );

% spirograph curve
[Time, WhCtrPos, MarkerPos, MarkerAngle] = ...
  AllBeziers( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

% time, parametrized by the arc of the marker or the ar of the wheel center
TimeFromMarker = zeros(1,size(MarkerPos,2));
TimeFromMarker(2:end) = cumsum( vecnorm( diff(MarkerPos,1,2), 2, 1 ) );

TimeFromWheel = zeros(1,size(WhCtrPos,2));
TimeFromWheel(2:end) = cumsum( vecnorm( diff(WhCtrPos,1,2), 2, 1 ) );

% patch
MarkerPos(:,end+1) = MarkerPos(:,1);
MarkerAngle(end+1) = MarkerAngle(1);

%%
% plot everything
figure()
plot(BezierPos(1,:),BezierPos(2,:),'blue')
hold on
%plot(WhCtrPos(1,:),WhCtrPos(2,:),'green')
axis equal
grid on
plot(MarkerPos(1,:),MarkerPos(2,:),'yellow')

%%
% make an animation
fps = 30;
MaxTime = ceil(TimeFromMarker(end)*fps)/fps;
nTimes = MaxTime*fps;

idxx = 1:size(TimeFromMarker,2);

%M(nTimes) = struct('cdata',[],'colormap',[]);

v = VideoWriter("test.mp4", 'MPEG-4');
open(v)

close all

f = figure();
f.Visible = 'off';

hold on
axis equal
axis off
plot(BezierPos(1,:),BezierPos(2,:),'blue')
ylim([-1.2,1.2])
%grid on

for i = 1:nTimes
  disp(i/nTimes)
  CurrPts = idxx(and(TimeFromMarker>=(i-1.1)/fps,TimeFromMarker<=(i+0.1)/fps));
  plot(MarkerPos(1,CurrPts),MarkerPos(2,CurrPts),'yellow')

  %M(i) = getframe;
  %frame = getframe(gcf);
  writeVideo(v,getframe)
end

hold off

f.Visible = 'on';

%movie(M)

close(v)

%%
% plot everything v2

% artistic choice
% video will last 30 seconds
TimeFromWheel = TimeFromWheel*(30/TimeFromWheel(end));

% parameters
fps = 30;
%MaxTime = ceil(TimeFromMarker(end)*fps)/fps;
MaxTime = ceil(TimeFromWheel(end)*fps)/fps;
nTimes = MaxTime*fps;

% help
%idxx = 1:size(TimeFromMarker,2);
idxx = 1:size(TimeFromWheel,2);
aang = 0:(2*pi/ ceil( 2*pi/(MaxDistDelta/WheelRadius) )):(2*pi);
circ = WheelRadius*[cos(aang); sin(aang)];
aux_angles = 0:(pi/3):(2*pi);

% original figure
f1 = figure(1);
hold on
axis equal
axis off
xlim([min(MarkerPos(1,:)) max(MarkerPos(1,:))])
ylim([min(MarkerPos(2,:)) max(MarkerPos(2,:))])
f1.Visible = 'off';
%
figure(1)
fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
%
f2 = figure(2);

% video object
v = VideoWriter("test_250908_11.mp4",'MPEG-4');
v.Quality = 100;
v.VideoBitsPerPixel
open(v)

% main loop
for i = 1:nTimes
  disp(i/nTimes)
  %CurrPts = idxx(and(TimeFromMarker>=(i-1.1)/fps,TimeFromMarker<=(i+0.1)/fps));
  CurrPts = idxx(and(TimeFromWheel>=(i-1.1)/fps,TimeFromWheel<=(i+0.1)/fps));
  %
  % add a few strokes of the marker, then copy to figure 2
  figure(1)
  plot(MarkerPos(1,CurrPts),MarkerPos(2,CurrPts),'red')
  clf(2)
  copyobj(gca,f2)
  %
  % add only the rotating wheel and the marker point
  figure(2)
  j = max(CurrPts);
  fill(WhCtrPos(1,j)+circ(1,:),WhCtrPos(2,j)+circ(2,:), 'b', 'EdgeColor', 'none','FaceAlpha',0.5); 
  for w = 1:size(aux_angles,2)
    plot(WhCtrPos(1,j)'+[0,cos(aux_angles(w)+MarkerAngle(j))*WheelRadius],WhCtrPos(2,j)'+[0,sin(aux_angles(w)+MarkerAngle(j))*WheelRadius],...
      'Color',[0,0,0,0.5])
  end
  scatter(MarkerPos(1,j),MarkerPos(2,j),10,'red','filled')

  writeVideo(v,getframe)
end

% finalize the video object and close the figure
close(v);
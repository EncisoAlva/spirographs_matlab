% the code for the curve is known to work, this should be the first real
% example
% the sketch of the moon took me some time, as I did it analytically

% designer stuff
MarkerAngle0 = pi;

% technical stuff
MaxDistDelta = 0.001;
CloseTol = 0.001;
MaxSpins = 20;

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
WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/9;


WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/15;
MarkerRadius = WheelRadius*(1);

% Bezier curve
BezierPos = AllBezierEval( CtrlPtsArray, MaxDistDelta );

% spirograph curve
[Time, WhCtrPos, MarkerPos, MarkerAngle] = ...
  AllBeziers( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

TimeFromMarker = zeros(1,size(MarkerPos,2));
TimeFromMarker(2:end) = cumsum( vecnorm( diff(MarkerPos,1,2), 2, 1 ) );

% patch
MarkerPos(:,end+1) = MarkerPos(:,1);
MarkerAngle(end+1) = MarkerAngle(1);

% plot everything
figure()
plot(BezierPos(1,:),BezierPos(2,:),'blue')
hold on
%plot(WhCtrPos(1,:),WhCtrPos(2,:),'green')
axis equal
grid on
plot(MarkerPos(1,:),MarkerPos(2,:),'yellow')

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
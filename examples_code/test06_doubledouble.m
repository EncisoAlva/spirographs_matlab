% the code for the curve is known to work, this should be the first real
% example
% the sketch of the moon took me some time, as I did it analytically

%%
% Bezier curves used, I made a small collection so far

% semicircles S
CtrlPtsArray = {[...
  [0,-1]',...
  [-(8/3)*tan(pi/4)',-1]',...
  [-(8/3)*tan(pi/4)',1]',...
  [0,1]'...
  ],[...
  [0,1]',...
  [-(4/3)*tan(pi/4)',1]',...
  [-(4/3)*tan(pi/4)',0]',...
  [0,0]'...
  ],[...
  [0,0]',...
  [(8/3)*tan(pi/4)',0]',...
  [(8/3)*tan(pi/4)',-2]',...
  [0,-2]'...
  ],[...
  [0,-2]'...
  [(4/3)*tan(pi/4)',-2]',...
  [(4/3)*tan(pi/4)',-1]',...
  [0,-1]'...
  ]};

% heart

CtrlPtsArray = {[...
  [-2,0]',...
  [-2,(4/3)*tan(pi/4)]',...
  [0,(4/3)*tan(pi/4)]',...
  [0,0]'...
  ],[...
  [0,0]',...
  [0,(4/3)*tan(pi/4)]',...
  [2,(4/3)*tan(pi/4)]',...
  [2,0]',...
  ],[...
  [2,0]',...
  [2,-2*(4/3)*tan(pi/16)]',...
  [2/sqrt(2),-2/sqrt(2)]'+2*[1/sqrt(2),1/sqrt(2)]'*(4/3)*tan(pi/16),...
  [2/sqrt(2),-2/sqrt(2)]'...
  ],[...
  [2/sqrt(2),-2/sqrt(2)]',...
  [2/sqrt(2),-2/sqrt(2)]'+[-.5,-.5]',...
  [0,-4/sqrt(2)]'+[.5,.5]',...
  [0,-4/sqrt(2)]'...
  ],[...
  [0,-4/sqrt(2)]',...
  [0,-4/sqrt(2)]'+[-.5,.5]',...
  [-2/sqrt(2),-2/sqrt(2)]'+[.5,-.5]',...
  [-2/sqrt(2),-2/sqrt(2)]',...
  ],[...
  [-2/sqrt(2),-2/sqrt(2)]'...
  [-2/sqrt(2),-2/sqrt(2)]'+2*[-1/sqrt(2),1/sqrt(2)]'*(4/3)*tan(pi/16),...
  [-2,-2*(4/3)*tan(pi/16)]',...
  [-2,0]'...
  ]};

% simple 4-spike star
[c1, c2] = HalfBezierSingle([...
  [-1,0]',...
  [-1,0]'+[1, 0]'*(4/3)*tan(pi/8), ...
  [ 0,1]'+[0,-1]'*(4/3)*tan(pi/8), ...
  [0,1]'
  ]);
CtrlPtsArray = {c2,[...
  [0,1]',...
  [0,1]'+[0,-1]'*(4/3)*tan(pi/8), ...
  [1,0]'+[-1,0]'*(4/3)*tan(pi/8), ...
  [1,0]'
  ],[...
  [1,0]',...
  [1,0]'+[-1,0]'*(4/3)*tan(pi/8), ...
  [0,-1]'+[0,1]'*(4/3)*tan(pi/8), ...
  [0,-1]'
  ],[...
  [0,-1]',...
  [0,-1]'+[0,1]'*(4/3)*tan(pi/8), ...
  [-1,0]'+[1,0]'*(4/3)*tan(pi/8), ...
  [-1,0]'
  ], c1};

% rotated square
CtrlPtsArray = {...
  LineToBezier([-1,1]',[0,2]'), ...
  LineToBezier([0,2]',[2,0]'), ...
  LineToBezier([2,0]',[0,-2]'), ...
  LineToBezier([0,-2]',[-2,0]'), ...
  LineToBezier([-2,0]',[-1,1]') ...
  };

%%
% self-intersecting shcape
xx = 1;
yy = 3;
[c1, c2] = HalfBezierSingle([1,0; 1+1+xx, -yy; -1-1-xx, -yy; -1,0]');
CtrlPtsArray = {c2,[...
  -1,0; 1+2*xx, 2*yy; -1-2*xx, 2*yy; 1, 0 ...
  ]',...
  c1};

if false
  CtrlPtsArray = FlipBezierAll(CtrlPtsArray);
end

% self-intersecting circle
CtrlPtsArray = {[...
  [0,-1]',...
  [0,-1]'+[-1,0]'*(4/3)*tan(pi/8),...
  [-1,0]'+[0,-1]'*(4/3)*tan(pi/8),...
  [-1,0]'...
  ],[...
  [-1,0]',...
  [-1,0]'+[0,1]'*(4/3)*tan(pi/8),...
  [1,2]'+[0,-1]'*(4/3)*tan(pi/8),...
  [1,2]'
  ],[...
  [1,2]',...
  [1,2]'+[0,1]'*(4/3)*tan(pi/8),...
  [0,3]'+[1,0]'*(4/3)*tan(pi/8),...
  [0,3]'
  ],[...
  [0,3]',...
  [0,3]'+[-1,0]'*(4/3)*tan(pi/8),...
  [-1,2]'+[0,1]'*(4/3)*tan(pi/8),...
  [-1,2]'...
  ],[...
  [-1,2]',...
  [-1,2]'+[0,-1]'*(4/3)*tan(pi/8),...
  [1,0]'+[0,1]'*(4/3)*tan(pi/8),...
  [1,0]'...
  ],[...
  [1,0]',...
  [1, 0]'+[0,-1]'*(4/3)*tan(pi/8),...
  [0,-1]'+[1, 0]'*(4/3)*tan(pi/8),...
  [0,-1]'
  ]};
CtrlPtsArray = FlipBezierAll(CtrlPtsArray);

% up V with circle borders
R = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];
CtrlPtsArray = {...
  R*LineToBezier([0,1]', [1,1]'), ...
  R*LineToBezier([1,1]', [1,-1]'), ...
  R*[...
  [1,-1]',...
  [1,-1]'+[0,-1]'*(2/3)*tan(pi/4),...
  [0,-1]'+[0,-1]'*(2/3)*tan(pi/4),...
  [0,-1]'
  ],...
  R*LineToBezier([0,-1]', [0,0]'), ...
  R*LineToBezier([0,0]', [-1,0]'), ...
  R*[...
  [-1,0]',...
  [-1,0]'+[-1,0]'*(2/3)*tan(pi/4),...
  [-1,1]'+[-1,0]'*(2/3)*tan(pi/4),...
  [-1,1]'
  ],...
  R*LineToBezier([-1,1]', [0,1]'), ...
  };

%%
%star with 5 spikes
aang   = -(0:(2*pi/5):(2*pi));
aang(end) = [];
PtsOut = [cos(aang + pi/2); sin(aang + pi/2)]*2;
PtsInn = [cos(aang + pi/2 - pi/5); sin(aang + pi/2 - pi/5)];

PtsAll = zeros(2,2*size(PtsInn,2)+1);
for i = 1:size(PtsInn,2)
  PtsAll(:,2*i-1) = PtsOut(:,i);
  PtsAll(:,2*i  ) = PtsInn(:,i);
end
PtsAll(:,end) = PtsOut(:,1);

CtrlPtsArray = {};
for i = 1:(size(PtsAll,2)-1)
  CurrCurve = zeros(2,4);
  CurrCurve(:,1) = PtsAll(:,i);
  CurrCurve(:,2) = (2/3)*PtsAll(:,i)+(1/3)*PtsAll(:,i+1);
  CurrCurve(:,3) = (1/3)*PtsAll(:,i)+(2/3)*PtsAll(:,i+1);
  CurrCurve(:,4) = PtsAll(:,i+1);
  %
  CtrlPtsArray{end+1} = CurrCurve;
end

CtrlPtsArray_backup = CtrlPtsArray;
[mid1, mid2] = HalfBezierSingle( CtrlPtsArray{1} );
CtrlPtsArray = { mid1, mid2 };
for i = 2:size(CtrlPtsArray_backup,2)
  CtrlPtsArray{end+1} = CtrlPtsArray_backup{i};
end

%%
figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray,2)
scatter(CtrlPtsArray{i}(1,:), CtrlPtsArray{i}(2,:))
end

BezOG  = AllBezierEval(CtrlPtsArray, 0.001);
figure()
hold on
axis equal
grid on
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');

%%
% parameters

% technical stuff
MaxDistDelta = 0.001;
CloseTol = 0.01;
MaxSpins = 100;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 5+5/7;
WheelMarkerRatio = 1;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.005*BezierArea(CtrlPtsArray, MaxDistDelta)/(pi));

%% 
% remove inner corners
if false
[CtrlPtsArray_rounded_flipped] = ...
  RemoveAllCorners( FlipBezierAll(CtrlPtsArray), CornerRoundingRadius, MaxDistDelta, false );
CtrlPtsArray = FlipBezierAll(CtrlPtsArray_rounded_flipped);

figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray,2)
scatter(CtrlPtsArray{i}(1,:), CtrlPtsArray{i}(2,:))
end
end

%% 
% remove corners inside and outside

WheelRadiusTol = 0.000001;

WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_tmp_inv] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray), WheelRadius_new, MaxDistDelta, true );
  [CtrlPtsArray_tmp] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray_tmp_inv), WheelRadius_new, MaxDistDelta, false );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_tmp,0.00001)/(2*pi))/WheelBezRatio
end
%CtrlPtsArray = CtrlPtsArray_tmp;
WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;

%CtrlPtsArray_backup = CtrlPtsArray;

BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);

figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray_tmp,2)
scatter(CtrlPtsArray_tmp{i}(1,:), CtrlPtsArray_tmp{i}(2,:))
end

%%
% checl
if false
BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);
BezNew = AllBezierEval(CtrlPtsArray_tmp, MaxDistDelta);

figure()
hold on
axis equal
grid on
fill(BezOG(1,:),BezOG(2,:), 'y', 'EdgeColor', 'none');
fill(BezNew(1,:),BezNew(2,:), 'r', 'EdgeColor', 'none');


[BezierPos, ~, ...
  ~,...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  ~,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray_tmp, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

figure()
fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')

BezOG = AllBezierEval(CtrlPtsArray, MaxDistDelta);
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none'); 

end

%%
[BezierPos1A, BezierPos2A, ...
  LocTime1A,...
  WhCtrPos1A, MarkerPos1A, MarkerAngle1A,...
  LocTime2A,...
  WhCtrPos2A, MarkerPos2A, MarkerAngle2A] = ...
  SetupCurves_2pts( CtrlPtsArray_tmp, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

[BezierPos1B, BezierPos2B, ...
  LocTime1B,...
  WhCtrPos1B, MarkerPos1B, MarkerAngle1B,...
  LocTime2B,...
  WhCtrPos2B, MarkerPos2B, MarkerAngle2B] = ...
  SetupCurves_2pts( FlipBezierAll(CtrlPtsArray_tmp), WheelRadius, MarkerRadius, MarkerAngle0+pi, ...
    MaxDistDelta, CloseTol, MaxSpins);

BezierPos1B   = flip(BezierPos1B, 2);
WhCtrPos1B    = flip( WhCtrPos1B, 2 );
MarkerPos1B   = flip(MarkerPos1B, 2);
MarkerAngle1B = flip(MarkerAngle1B, 2);
%
BezierPos2B   = flip(BezierPos2B, 2);
WhCtrPos2B    = flip(WhCtrPos2B, 2);
MarkerPos2B   = flip(MarkerPos2B, 2);
MarkerAngle2B = flip(MarkerAngle2B, 2);

LocTime1B = max(LocTime1B) - flip(LocTime1B);
LocTime2B = max(LocTime2B) - flip(LocTime2B);

%%
% preview
figure()
fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(MarkerPos1A(1,:),MarkerPos1A(2,:),'yellow')
plot(MarkerPos2A(1,:),MarkerPos2A(2,:),'magenta')
plot(MarkerPos1B(1,:),MarkerPos1B(2,:),'yellow')
plot(MarkerPos2B(1,:),MarkerPos2B(2,:),'magenta')

%%
% video parameters
TotalTime = 60;
AfterTime = 5;

VidName = 'doubledouble250916_2';

%%
% video

close all

TimeFromBez1A = zeros(1,size(BezierPos1A,2));
TimeFromBez1A(2:end) = cumsum( vecnorm( diff(BezierPos1A,1,2), 2, 1 ) );

TimeFromBez2A = zeros(1,size(BezierPos2A,2));
TimeFromBez2A(2:end) = cumsum( vecnorm( diff(BezierPos2A,1,2), 2, 1 ) );

TimeFromBez1B = zeros(1,size(BezierPos1B,2));
TimeFromBez1B(2:end) = cumsum( vecnorm( diff(BezierPos1B,1,2), 2, 1 ) );

TimeFromBez2B = zeros(1,size(BezierPos2B,2));
TimeFromBez2B(2:end) = cumsum( vecnorm( diff(BezierPos2B,1,2), 2, 1 ) );

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
f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off
xlim([min(MarkerPos1A(1,:)) max(MarkerPos1A(1,:))])
ylim([min(MarkerPos1A(2,:)) max(MarkerPos1A(2,:))])
xlim([min(MarkerPos2A(1,:)) max(MarkerPos2A(1,:))])
ylim([min(MarkerPos2A(2,:)) max(MarkerPos2A(2,:))])
%
%fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
plot(BezNew(1,:),BezNew(2,:),'k','LineWidth',2)
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
  plot(MarkerPos1A(1,CurrPts1A),MarkerPos1A(2,CurrPts1A),'magenta')
  plot(MarkerPos2A(1,CurrPts2A),MarkerPos2A(2,CurrPts2A),'yellow')
  plot(MarkerPos1B(1,CurrPts1B),MarkerPos1B(2,CurrPts1B),'yellow')
  plot(MarkerPos2B(1,CurrPts2B),MarkerPos2B(2,CurrPts2B),'magenta')
  %
  j1A = max(CurrPts1A);
  j2A = max(CurrPts2A);
  j1B = max(CurrPts1B);
  j2B = max(CurrPts2B);
  if ~isempty(j1A)
    RefWheelCtrA = WhCtrPos1A(:,j1A);
    RefAngleA = MarkerAngle1A(j1A);
  else
    RefWheelCtrA = WhCtrPos2A(:,j2A);
    RefAngleA = MarkerAngle2A(j2A);
  end
  if ~isempty(j1B)
    RefWheelCtrB = WhCtrPos1B(:,j1B);
    RefAngleB = MarkerAngle1B(j1B);
  else
    RefWheelCtrB = WhCtrPos2B(:,j2B);
    RefAngleB = MarkerAngle2B(j2B);
  end
  fill(RefWheelCtrA(1)+circ(1,:),RefWheelCtrA(2)+circ(2,:), 'b', 'EdgeColor', 'none','FaceAlpha',0.5); 
  fill(RefWheelCtrB(1)+circ(1,:),RefWheelCtrB(2)+circ(2,:), 'b', 'EdgeColor', 'none','FaceAlpha',0.5); 
  for w = 1:size(aux_angles,2)
    plot(RefWheelCtrA(1)+[0,cos(aux_angles(w)+RefAngleA)*WheelRadius],RefWheelCtrA(2)+[0,sin(aux_angles(w)+RefAngleA)*WheelRadius],...
      'Color',[0,0,0,0.5])
    plot(RefWheelCtrB(1)+[0,cos(aux_angles(w)+RefAngleB)*WheelRadius],RefWheelCtrB(2)+[0,sin(aux_angles(w)+RefAngleB)*WheelRadius],...
      'Color',[0,0,0,0.5])
  end
  scatter(MarkerPos1A(1,j1A),MarkerPos1A(2,j1A),10,'magenta','filled')
  scatter(MarkerPos2A(1,j2A),MarkerPos2A(2,j2A),10,'yellow','filled')
  scatter(MarkerPos1B(1,j1B),MarkerPos1B(2,j1B),10,'yellow','filled')
  scatter(MarkerPos2B(1,j2B),MarkerPos2B(2,j2B),10,'magenta','filled')
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
% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','LetterC'));
CtrlPtsArray = CtrlPtsArray{1};

%%
% load from file
AllCtrlPtsArray = LoadSVG( './curves_svg/CopperBlack_M.svg' );
CtrlPtsArray = AllCtrlPtsArray{1};

%%
% pre-processing
CtrlPtsArray = RemovePointCurves( CtrlPtsArray, 0.0001 );

% this one is because I'm using absolute tolerance instead of relative
CtrlPtsArray = RescaleShape( CtrlPtsArray, 2, 2 );

% line with bad encoding, the normal vector will be wrong
nCurves = size(CtrlPtsArray,2);
for i = 1:nCurves
  CurrCurve = CtrlPtsArray{i};
  if norm(CurrCurve(:,1)-CurrCurve(:,2)) < 0.001
    if abs( norm(CurrCurve(:,2)-CurrCurve(:,3)) + norm(CurrCurve(:,3)-CurrCurve(:,4)) - norm(CurrCurve(:,2)-CurrCurve(:,4)) ) < 0.001
      CtrlPtsArray{i} = LineToBezier( CurrCurve(:,1), CurrCurve(:,4) );
      continue
    end
  end
  if norm(CurrCurve(:,3)-CurrCurve(:,4)) < 0.001
    if abs( norm(CurrCurve(:,1)-CurrCurve(:,2)) + norm(CurrCurve(:,2)-CurrCurve(:,3)) - norm(CurrCurve(:,1)-CurrCurve(:,3)) ) < 0.001
      CtrlPtsArray{i} = LineToBezier( CurrCurve(:,1), CurrCurve(:,4) );
      continue
    end
  end
end

% rotate half a spin
for i = 1:size(CtrlPtsArray, 2)
  CtrlPtsArray{i} = [1,0; 0,-1] * CtrlPtsArray{i};
end

% rotate by an angle
if false
  th = pi/2;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  for i = 1:size(CtrlPtsArray, 2)
    CtrlPtsArray{i} = ROT * CtrlPtsArray{i};
  end
end

% invert orientation, if needed
if false
  CtrlPtsArray = FlipBezierAll(CtrlPtsArray);
end

%%
% show control points
figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray,2)
  scatter(CtrlPtsArray{i}(1,:), CtrlPtsArray{i}(2,:))
end

% show shape
BezOG  = AllBezierEval(CtrlPtsArray, 0.001);
figure()
hold on
axis equal
grid on
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');

%%
% parameters

% technical stuff
MaxDistDelta = 0.005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.000001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 1;
WheelMarkerRatio = 4/5;

Shift  = 0;
Halfen = true;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.005*BezierArea(CtrlPtsArray, MaxDistDelta)/(pi));

%% 
% remove inner corners
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

%% 
% remove outer corners
WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio
while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_new] = ...
    RemoveAllCorners( CtrlPtsArray, WheelRadius_new, MaxDistDelta, true );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_new,0.00001)/(2*pi))/WheelBezRatio
end
%CtrlPtsArray = CtrlPtsArray_new;
WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;
CtrlPtsArray_new = ShiftBezierAll( CtrlPtsArray_new, Shift, Halfen );

if false
  CtrlPtsArray_new = CtrlPtsArray;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio
  WheelRadius  = WheelRadius_new;
  MarkerRadius = WheelRadius*WheelMarkerRatio;
end

%% 

% show control points
figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray_new,2)
scatter(CtrlPtsArray_new{i}(1,:), CtrlPtsArray_new{i}(2,:))
end

% difference from rounding
BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);
figure()
hold on
axis equal
grid on
fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none');
%fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');

%%
ColorVector = {'yellow','magenta', 'red', 'red'};

%%

% preview result
[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);

figure()
hold on
axis equal
grid on
fill(BezNew(1,:),BezNew(2,:), .15*[1,1,1], 'EdgeColor', 'none')
plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),ColorVector{1}, 'LineWidth',2)
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),ColorVector{2}, 'LineWidth',2)
set(gca,'color', 'k');

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 +pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),ColorVector{3}, 'LineWidth',2)
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),ColorVector{4}, 'LineWidth',2)
set(gca,'color', 'k');

%%
[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);

figure()
hold on
axis equal
axis off
grid off
plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),ColorVector{1}, 'LineWidth',2)
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),ColorVector{2}, 'LineWidth',2)
set(gca,'color', 'k');

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 +pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),ColorVector{3}, 'LineWidth',2)
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),ColorVector{4}, 'LineWidth',2)
set(gca,'color', 'k');

%%
figure()
hold on
axis equal
axis off
grid off
fill(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:), 'magenta', 'EdgeColor', 'none');


%%

% make curves
[ DecorativeBez,...
  AllBezierPos1, AllLocTime1, ...
  AllWhCtrPos1, AllMarkerPos1, AllMarkerAngle1 ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);
[ ~,...
  AllBezierPos2, AllLocTime2, ...
  AllWhCtrPos2, AllMarkerPos2, AllMarkerAngle2 ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0+pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

%
AllBezierPos   = [ AllBezierPos1,   AllBezierPos2 ];
AllWhCtrPos    = [ AllWhCtrPos1,    AllWhCtrPos2  ];
AllMarkerPos   = [ AllMarkerPos1,   AllMarkerPos2 ];
AllMarkerAngle = [ AllMarkerAngle1, AllMarkerAngle2 ];
AllLocTime     = [ AllLocTime1,     AllLocTime2 ];

% video
MakeVideo_4pts( WheelRadius, 'Wheel', 'in', ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector,...
  MaxDistDelta, ...
  40, 10, 'test_251002_01' )

%  {'yellow', 'magenta', 'red', 'red'},...

%%

% make curves
[ DecorativeBez,...
  AllBezierPos, ~, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

% video
MakeVideo_2pts( WheelRadius, 'Bezier', 'in', ...
  DecorativeBez,...
  AllBezierPos, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  MaxDistDelta, ...
  40, 10, 'test_251001_02' )

%Wheel
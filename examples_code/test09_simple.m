% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','Astroid4'));
CtrlPtsArray = CtrlPtsArray{1};

%%
% load from file
AllCtrlPtsArray = LoadSVG( './curves_svg/coin_line.svg' );
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
  th = -pi/2;
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

PlotBezierCtrlPts(CtrlPtsArray)

%CtrlPtsArray = ShiftBezierAll( CtrlPtsArray, 15, false );
%CtrlPtsArray_back = CtrlPtsArray;
%CtrlPtsArray = CtrlPtsArray_back;

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

WheelBezRatio = 4;
WheelMarkerRatio = 4/5;

Shift  = 0;
Halfen = false;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.0001*BezierArea(CtrlPtsArray, MaxDistDelta)/(pi));

%% 
% remove inner corners
[CtrlPtsArray_rounded_flipped] = ...
  RemoveAllCorners( FlipBezierAll(CtrlPtsArray), CornerRoundingRadius, MaxDistDelta, false );
CtrlPtsArray_tmp = FlipBezierAll(CtrlPtsArray_rounded_flipped);

PlotBezierCtrlPts(CtrlPtsArray_tmp)

%
CtrlPtsArray = CtrlPtsArray_tmp;

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
PlotBezierCtrlPts(CtrlPtsArray_new)

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

%ColorVector = {'yellow','white', 'red', 'red'};

%ColorVector = {'red','white', 'red'};

%CtrlPtsArray_new = ShiftBezierAll( CtrlPtsArray_new, -3, false );

%%

CurveOpts = {};
CurveOpts.CloseEnds = false;

% preview result
[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, CurveOpts);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);

figure()
hold on
axis equal
grid on
fill(BezNew(1,:),BezNew(2,:), .15*[1,1,1], 'EdgeColor', 'none')
plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),ColorVector{1}, 'LineWidth',2)
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),ColorVector{2}, 'LineWidth',2)
set(gca,'color', 'k');
scatter(AllMarkerPos{1}(1,1),AllMarkerPos{1}(2,1),'red','filled','o')
scatter(AllMarkerPos{2}(1,1),AllMarkerPos{2}(2,1),'red','filled','o')

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 +pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, CurveOpts);

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

XL = xlim;
YL = ylim;
%
x_ran = XL(2) - XL(1);
y_ran = YL(2) - YL(1);
ref_ran = max(x_ran, y_ran);
%
XL(1) = XL(1) - (ref_ran - x_ran)/2;
XL(2) = XL(2) + (ref_ran - x_ran)/2;
YL(1) = YL(1) - (ref_ran - y_ran)/2;
YL(2) = YL(2) + (ref_ran - y_ran)/2;
%
xlim([XL(1) XL(2)])
ylim([YL(1) YL(2)])

%%
figure()
hold on
axis equal
axis off
grid off

fill(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:), 'magenta', 'EdgeColor', 'none');
fill(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:), 'magenta', 'EdgeColor', 'none');

%%

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 +pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

figure()
hold on
axis equal
axis off
grid off
AlternatePos = [AllMarkerPos{1}, AllMarkerPos{2}];

fill(AlternatePos(1,:),AlternatePos(2,:), 'red', 'EdgeColor', 'none');


%%

figure()
hold on
axis equal
axis off
grid off

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 +pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);
AlternatePos = [AllMarkerPos{1}, AllMarkerPos{2}];
fill(AlternatePos(1,:),AlternatePos(2,:), 'green', 'EdgeColor', 'none');

[ ~, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0 , ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);
AlternatePos = [AllMarkerPos{1}, AllMarkerPos{2}];
fill(AlternatePos(1,:),AlternatePos(2,:), 'red', 'EdgeColor', 'none');

XL = xlim;
YL = ylim;
%
x_ran = XL(2) - XL(1);
y_ran = YL(2) - YL(1);
ref_ran = max(x_ran, y_ran);
%
XL(1) = XL(1) - (ref_ran - x_ran)/2;
XL(2) = XL(2) + (ref_ran - x_ran)/2;
YL(1) = YL(1) - (ref_ran - y_ran)/2;
YL(2) = YL(2) + (ref_ran - y_ran)/2;
%
xlim([XL(1) XL(2)])
ylim([YL(1) YL(2)])


%%
CurveOpts = {};
CurveOpts.CloseEnds = false;

% make curves
[ DecorativeBez,...
  AllBezierPos1, AllLocTime1, ...
  AllWhCtrPos1, AllMarkerPos1, AllMarkerAngle1 ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, CurveOpts);
[ ~,...
  AllBezierPos2, AllLocTime2, ...
  AllWhCtrPos2, AllMarkerPos2, AllMarkerAngle2 ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0+pi/2, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, CurveOpts);

%
AllBezierPos   = [ AllBezierPos1,   AllBezierPos2 ];
AllWhCtrPos    = [ AllWhCtrPos1,    AllWhCtrPos2  ];
AllMarkerPos   = [ AllMarkerPos1,   AllMarkerPos2 ];
AllMarkerAngle = [ AllMarkerAngle1, AllMarkerAngle2 ];
AllLocTime     = [ AllLocTime1,     AllLocTime2 ];

ExtraOpts = {};
ExtraOpts.Plot2Circles = false;
ExtraOpts.Format = 'mp4';
ExtraOpts.Orientation = 'in';
ExtraOpts.Ratio = 16/9;
ExtraOpts.TimerefCurve = 'Wheel';
ExtraOpts.LineWidth = 2;

% video
MakeVideo_4pts( WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector,...
  MaxDistDelta, ...
  40, 10, 'test_251010_01', ExtraOpts )

%  {'yellow', 'magenta', 'red', 'red'},...

%%

CurveOpts = {};
CurveOpts.CloseEnds = false;

% make curves
[ DecorativeBez,...
  AllBezierPos, ~, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, CurveOpts);

% extra options
ExtraOpts = {};
ExtraOpts.Plot2Circles = false;
ExtraOpts.Format = 'mp4';
ExtraOpts.Orientation = 'in';
ExtraOpts.Ratio = 16/9;
ExtraOpts.TimerefCurve = 'Wheel';
ExtraOpts.LineWidth = 2;

% video
MakeVideo_2pts( WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector,...
  MaxDistDelta, ...
  40, 10, 'test_251010_00', ExtraOpts )

%Wheel
%Bezier
%Average
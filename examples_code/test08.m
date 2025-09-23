% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','BumpCircle'));
CtrlPtsArray = CtrlPtsArray{1};

%%
% load from file
AllCtrlPtsArray = LoadSVG( './curves_svg/cinquefoil.svg' );
CtrlPtsArray = AllCtrlPtsArray{1};

%%
% pre-processing
CtrlPtsArray = RemovePointCurves( CtrlPtsArray, 0.0001 );

% this one is because I'm using absolute tolerance instead of relative
CtrlPtsArray = RescaleShape( CtrlPtsArray, 2, 2 );

%CtrlPtsArray = FlipBezierAll(CtrlPtsArray);

if false
for i = 1:size(CtrlPtsArray, 2)
  CtrlPtsArray{i} = [1,0; 0,-1] * CtrlPtsArray{i};
end
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

WheelBezRatio = 40+1/3;
WheelMarkerRatio = 4/5;

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
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;
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

%CtrlPtsArray_new = CtrlPtsArray;

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
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');

% preview result
[~, ~, ~, ~, MarkerPos1, ~, ~, ~, MarkerPos2, ~] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);

figure()
hold on
axis equal
grid on
fill(BezNew(1,:),BezNew(2,:), .15*[1,1,1], 'EdgeColor', 'none')
%fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none'); 
%fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none'); 
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')
set(gca,'color', 'k');

figure()
hold on
axis equal
grid off
axis off
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')
set(gca,'color', 'k');

figure()
hold on
axis equal
axis off
grid off
%fill(MarkerPos2(1,:),MarkerPos2(2,:), 'yellow', 'EdgeColor', 'none');
fill(MarkerPos1(1,:),MarkerPos1(2,:), 'magenta', 'EdgeColor', 'none');


%%

% make curves
[BezierPos, ~, ...
  ~,...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  ~,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

% video
close all
MakeVideo_2pts( WheelRadius, ...
  BezOG, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2,...
  MaxDistDelta,...
  60, 2, 5, 'test_250923_09' )
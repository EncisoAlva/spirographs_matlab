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

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 10+1/9;
WheelMarkerRatio = 4/5;

% specific to this example
%WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;
%MarkerRadius = WheelRadius*WheelMarkerRatio;

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
% remove corners

WheelRadiusTol = 0.000001;

WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_tmp] = ...
    RemoveAllCorners( CtrlPtsArray, WheelRadius_new, MaxDistDelta, true );
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
fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none');
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');


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
% video

close all
MakeVideo_2pts( WheelRadius, ...
  BezOG, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2,...
  MaxDistDelta,...
  60, 5, 'test_250918_15' )
%  60, 5, 'test_250914_02' )
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
WheelMarkerRatio = 1;

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

%%

VidName = 'moving_250923_23';

TotalTime = 10;
LoopTime0 = 3;


%%

% rescale
CtrlPtsArray_big = RescaleShape(CtrlPtsArray_new, 100, 100);
BezNew = AllBezierEval(CtrlPtsArray_big, MaxDistDelta);

%
fps = 30;

set(0, 'DefaultFigureColor', 'k');

% computed parameters
LoopTime  = TotalTime/ ceil(TotalTime/LoopTime0);
nLoops   = TotalTime/LoopTime;
AngleIncr = 2*pi/(fps*LoopTime);

ratio = flip([1/3, 2/3, 1, 4/3, 5/3, 6/3]);

close all

v = VideoWriter(strcat(VidName,".avi"),'Uncompressed AVI');
v.VideoBitsPerPixel
%v.Quality = 100;
open(v)

% main loop
WB = waitbar(0,strcat('Generating video (',VidName,'.mp4)...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2025)');

f1 = figure('Visible','off','Name','Just the curve');

TotFrames = nLoops*LoopTime*fps;
counter = 0;
for currLoop = 1:nLoops
for currFrame = 1:(fps*LoopTime)
waitbar(counter/TotFrames, WB, sprintf('Generating video (%.2f%%)', counter/TotFrames*100));

clf(f1)
hold on
axis equal
axis off
grid off
set(gca,'color', 'k');
xlim([min(BezNew(1,:))-WheelRadius*max(ratio) , max(BezNew(1,:))+WheelRadius*max(ratio)])
ylim([min(BezNew(2,:))-WheelRadius*max(ratio) , max(BezNew(2,:))+WheelRadius*max(ratio)])

% make curves for different marker radius ratios
for i = 1:size(ratio,2)
  [~, ~, ~, ~, MarkerPos1, ~, ~, ~, ~, ~] = ...
      SetupCurves_2pts(CtrlPtsArray_big, WheelRadius, MarkerRadius * ratio(i), MarkerAngle0+currFrame*AngleIncr, ...
      MaxDistDelta, CloseTol, MaxSpins);
  
  plot(MarkerPos1(1,:), MarkerPos1(2,:), 'Color', [1,0,1]*(i/size(ratio,2)));
end
writeVideo(v, getframe);

counter = counter+1;
end
end

close(v);
delete(WB);
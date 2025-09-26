% Two circles, one inside and one outside, marking a total of 4 points. The
% shape is rounded prior to accomodate for it.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','Heart'));
CtrlPtsArray = CtrlPtsArray{1};

%CtrlPtsArray = Fidget3;

%%
% load from file
AllCtrlPtsArray = LoadSVG( './curves_svg/wobbly01.svg' );
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

for i = 1:size(CtrlPtsArray, 2)
  CtrlPtsArray{i} = [1,0; 0,-1] * CtrlPtsArray{i};
end

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
MaxDistDelta = 0.0002;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.0001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 10 + 1/12;
WheelMarkerRatio = 1;

%% 
% remove corners inside and outside
WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_new_inv] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray), WheelRadius_new/3, MaxDistDelta, true );
  [CtrlPtsArray_new] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray_new_inv), WheelRadius_new, MaxDistDelta, false );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_new,0.00001)/(2*pi))/WheelBezRatio
end

WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;

if false

WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio
WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;

CtrlPtsArray_new = CtrlPtsArray;

end

%%

[ ~,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_4pts_smaller( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, 3);


%%
% visualize result

nCurves = size(CtrlPtsArray_new,2);
nSpins = max(AllLocTime{2})/nCurves;

figure()
hold on
axis equal
axis off
grid off

for i = 0:(nSpins-1)
  if ismember(i, [0,1, 4,5, 8,9])
    currColor = [1,0,1,.5];
  else
    currColor = [1,1,0,.5];
  end
  for j = 2:2
    currCurve = AllMarkerPos{j}(:,(floor(AllLocTime{j}/nCurves)>=i)&(floor(AllLocTime{j}/nCurves)<(i+1)));
    plot(currCurve(1,:),currCurve(2,:),'Color',currColor, 'LineWidth',0.2)
  end
  %
  %
  if ismember(i, [0,1])
    currColor = [1,0,1,.5];
  else
    currColor = [1,1,0,.5];
  end
   for j = 3:3
    currCurve = AllMarkerPos{j}(:,(floor(AllLocTime{j}/nCurves)>=i)&(floor(AllLocTime{j}/nCurves)<(i+1)));
    plot(currCurve(1,:),currCurve(2,:),'Color',currColor, 'LineWidth',0.2)
  end

end

%%

TotalTime = 5;
LoopTime0 = 5;

VidName = 'snake_250925_19';

%%

%
fps = 30;
BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);

% computed parameters
LoopTime  = TotalTime/ ceil(TotalTime/LoopTime0);
nLoops   = TotalTime/LoopTime;
AngleIncr = 2*pi/(fps*LoopTime);

close all
set(0, 'DefaultFigureColor', 'k');

v = VideoWriter(strcat(VidName,".mp4"),'MPEG-4');
v.Quality = 100;
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
set(f1, 'Color', 'k');
hold on
axis equal
axis off
grid off
xlim([min(BezNew(1,:))-2*WheelRadius , max(BezNew(1,:))+2*WheelRadius])
ylim([min(BezNew(2,:))-2*WheelRadius , max(BezNew(2,:))+2*WheelRadius])

[ ~,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_4pts_smaller( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0+currFrame*AngleIncr, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, 3);

% make curves for different marker radius ratios
for i = 0:(nSpins-1)
  if ismember(i, [0,1, 4,5, 8,9])
    currColor = [1,0,1,.5];
  else
    currColor = [1,1,0,.5];
  end
  for j = 2:2
    currCurve = AllMarkerPos{j}(:,(floor(AllLocTime{j}/nCurves)>=i)&(floor(AllLocTime{j}/nCurves)<(i+1)));
    plot(currCurve(1,:),currCurve(2,:),'Color',currColor, 'LineWidth',0.2)
  end
  %
  %
  if ismember(i, [0,1])
    currColor = [1,0,1,.5];
  else
    currColor = [1,1,0,.5];
  end
   for j = 3:3
    currCurve = AllMarkerPos{j}(:,(floor(AllLocTime{j}/nCurves)>=i)&(floor(AllLocTime{j}/nCurves)<(i+1)));
    plot(currCurve(1,:),currCurve(2,:),'Color',currColor, 'LineWidth',0.2)
  end

end


writeVideo(v, getframe);

counter = counter+1;
end
end

close(v);
delete(WB);
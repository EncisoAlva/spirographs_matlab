% Two circles, one inside and one outside, marking a total of 4 points. The
% shape is rounded prior to accomodate for it.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','BumpCircle'));
CtrlPtsArray = CtrlPtsArray{1};

%%
% 5 eyes thingy
aang = -(0:(2*pi/5):2*pi)+pi/2+pi/5;
hex_pts = [cos(aang); sin(aang)]*cot(pi/5);

CtrlPtsArray = {};
for i = 1:5
  CurrCtrl = zeros(2,4);
  CurrCtrl(:,1) = hex_pts(:,i  );
  CurrCtrl(:,2) = hex_pts(:,i  )*(1 +(4/3)*tan(7*pi/5) );
  CurrCtrl(:,3) = hex_pts(:,i+1)*(1 +(4/3)*tan(7*pi/5) );
  CurrCtrl(:,4) = hex_pts(:,i+1);
  CtrlPtsArray{end+1} = CurrCtrl;
end

% change starting point (artistic choice)
[c1,c2] = HalfBezierSingle(CtrlPtsArray{1});
CtrlPtsArray{1} = c2;
CtrlPtsArray{end+1} = c1;

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
MaxDistDelta = 0.0005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.0001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 15+1/5;
WheelMarkerRatio = 1;

%% 
% remove corners inside and outside
WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_new_inv] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray), WheelRadius_new, MaxDistDelta, true );
  [CtrlPtsArray_new] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray_new_inv), WheelRadius_new, MaxDistDelta, false );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_new,0.00001)/(2*pi))/WheelBezRatio
end
WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;

%%
% check if results are good enough

% control points[
figure()
hold on
axis equal
grid on
for i = 1:size(CtrlPtsArray_new,2)
scatter(CtrlPtsArray_new{i}(1,:), CtrlPtsArray_new{i}(2,:))
end

% difference induced by rounding
BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);
BezNew = AllBezierEval(CtrlPtsArray_new, MaxDistDelta);
figure()
hold on
axis equal
grid on
fill(BezOG(1,:), BezOG(2,:),  'y', 'EdgeColor', 'none');
fill(BezNew(1,:),BezNew(2,:), 'r', 'EdgeColor', 'none');

% preview of the result
[~, ~, ...
  ~,...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  ~,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);
figure()
fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')


%%
% preview
[ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_4pts( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

figure()
fill(DecorativeBez(1,:),DecorativeBez(2,:), 'k', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),'yellow')
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),'magenta')
plot(AllMarkerPos{3}(1,:),AllMarkerPos{3}(2,:),'yellow')
plot(AllMarkerPos{4}(1,:),AllMarkerPos{4}(2,:),'magenta')

%%
% video parameters
TotalTime = 60;
AfterTime = 5;

VidName = 'doubledouble250920_12';

%%
% update curve 
CtrlPtsArray = CtrlPtsArray_new;

% make curves
[ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_4pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

% video
MakeVideo_4pts( WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  MaxDistDelta, ...
  TotalTime, AfterTime, VidName )
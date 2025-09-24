% Two circles, one inside and one outside, marking a total of 4 points. The
% shape is rounded prior to accomodate for it.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
CtrlPtsArray = struct2cell(load('ExampleCurves.mat','LetterI'));
CtrlPtsArray = CtrlPtsArray{3};

%CtrlPtsArray = Fidget3;

%%
% load from file
AllCtrlPtsArray = LoadSVG( './curves_svg/MYNAME.svg' );
CtrlPtsArray = AllCtrlPtsArray{2};

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

WheelBezRatio = 20 + 1/6;
WheelMarkerRatio = 1;

%% 
% remove corners inside and outside
WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_new_inv] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray), WheelRadius_new/2, MaxDistDelta, true );
  [CtrlPtsArray_new] = ...
    RemoveAllCorners( FlipBezierAll(CtrlPtsArray_new_inv), WheelRadius_new, MaxDistDelta, false );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_new,0.00001)/(2*pi))/WheelBezRatio
end

WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;

%CtrlPtsArray_new = CtrlPtsArray;

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

%%
% preview
[ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
  SetupCurves_4pts_smaller( CtrlPtsArray_new, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta/2, CloseTol, MaxSpins,2);

figure()
%fill(DecorativeBez(1,:),DecorativeBez(2,:), [.4,.4,.4], 'EdgeColor', 'none'); 
hold on
axis equal
axis off
grid off
set(gcf,'Color','white')
%plot(AllMarkerPos{1}(1,:),AllMarkerPos{1}(2,:),'yellow')
plot(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:),'blue')
plot(AllMarkerPos{3}(1,:),AllMarkerPos{3}(2,:),'blue')
%plot(AllMarkerPos{4}(1,:),AllMarkerPos{4}(2,:),'yellow')

%%

figure()
%fill(DecorativeBez(1,:),DecorativeBez(2,:), [.4,.4,.4], 'EdgeColor', 'none'); 
hold on
axis equal
axis off
grid off
%fill(MarkerPos2(1,:),MarkerPos2(2,:), 'yellow', 'EdgeColor', 'none');
fill(AllMarkerPos{2}(1,:),AllMarkerPos{2}(2,:), 'magenta', 'EdgeColor', 'none');
fill(AllMarkerPos{3}(1,:),AllMarkerPos{3}(2,:), 'k', 'EdgeColor', 'none');

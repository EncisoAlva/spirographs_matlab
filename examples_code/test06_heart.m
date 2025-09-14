% the code for the curve is known to work, this should be the first real
% example
% the sketch of the moon took me some time, as I did it analytically

%%
% Bezier curves used, I made a small collection so far

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

if false
  CtrlPtsArray = FlipBezierAll(CtrlPtsArray);
end

%%
% parameters

% technical stuff
MaxDistDelta = 0.001;
CloseTol = 0.001;
MaxSpins = 100;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 33/7;
WheelMarkerRatio = 1;

% specific to this example
WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;
MarkerRadius = WheelRadius*WheelMarkerRatio;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.001*BezierArea(CtrlPtsArray, MaxDistDelta)/(pi));

%% do curves
if false
[BezierPos, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

% preview
figure()
fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')
end

%% 
% remove inner corners
[CtrlPtsArray_rounded_flipped] = ...
  RemoveAllCorners( FlipBezierAll(CtrlPtsArray), CornerRoundingRadius, MaxDistDelta, true );
CtrlPtsArray = FlipBezierAll(CtrlPtsArray_rounded_flipped);

%WheelRadius = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;

%[BezierPos, ...
%  WhCtrPos1, MarkerPos1, MarkerAngle1,...
%  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
%  SetupCurves_2pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
%    MaxDistDelta, CloseTol, MaxSpins);

% preview
%figure()
%fill(BezierPos(1,:),BezierPos(2,:), 'k', 'EdgeColor', 'none'); 
%hold on
%axis equal
%grid on
%plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
%plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')

%% 
% remove corners

WheelRadiusTol = 0.0000001;

WheelRadius_old = Inf;
WheelRadius_new = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;

while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [CtrlPtsArray_tmp] = ...
    RemoveAllCorners( CtrlPtsArray, WheelRadius_new, MaxDistDelta, false );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (BezierPerimeter(CtrlPtsArray_tmp,0.00001)/(2*pi))/WheelBezRatio;
end
%CtrlPtsArray = CtrlPtsArray_tmp;
WheelRadius = WheelRadius_new;

%CtrlPtsArray_backup = CtrlPtsArray;

%%
% checl
BezOG  = AllBezierEval(CtrlPtsArray, MaxDistDelta);
BezNew = AllBezierEval(CtrlPtsArray_tmp, MaxDistDelta);

figure()
hold on
axis equal
grid on
fill(BezNew(1,:),BezNew(2,:), 'y', 'EdgeColor', 'none');
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none');


[BezierPos, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray_tmp, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

figure()
fill(BezierPos(1,:),BezierPos(2,:), 'y', 'EdgeColor', 'none'); 
hold on
axis equal
grid on
plot(MarkerPos1(1,:),MarkerPos1(2,:),'yellow')
plot(MarkerPos2(1,:),MarkerPos2(2,:),'magenta')

BezOG = AllBezierEval(CtrlPtsArray, MaxDistDelta);
fill(BezOG(1,:),BezOG(2,:), 'r', 'EdgeColor', 'none'); 

%%
% video

close all
MakeVideo_2pts( WheelRadius, ...
  BezOG, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2,...
  MaxDistDelta,...
  60, 5, 'test_250914_02' )
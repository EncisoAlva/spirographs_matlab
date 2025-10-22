% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','Number8'));
BPath = BPath_pack{1};
clear BPath_pack

%%
% check available curves in the example file
who -file ExampleCollections.mat

% load curve
BPath_pack1 = struct2cell(load('ExampleCollections.mat','Circlegon'));
BPath_pack2 = BPath_pack1{1};
BPath = BPath_pack2{4};
clear BPath_pack1 BPath_pack2

%%
% load from file
BPatch_pack = LoadSVG( './curves_svg/coin_line.svg' );
BPath = BPath_pack{1};
clear BPath_pack

%%
% pre-processing
BPath = RemovePointCurves( BPath, 0.0001 );

% this one is because I'm using absolute tolerance instead of relative
BPath = RescalePath( BPath, 2, 2 );

% line with bad encoding, the normal vector will be wrong
BPath = ForceCubicLines( BPath );

% rotate half a spin
for i = 1:size(BPath, 2)
  BPath{i} = [1,0; 0,-1] * BPath{i};
end

% rotate by an angle
BPath = RotatePath( BPath, -pi/6 );

%%
% show control points

PlotPath(BPath)

BPath = ShiftPath( BPath, 2, false );

%%
% parameters

% technical stuff
Tol = 0.005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.000001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 5/2;
WheelMarkerRatio = 4/5;

%Shift  = -2;
%Halfen = true;
%Shift  = 3;
%Halfen = false;
Shift  = 0;
Halfen = false;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.001*PathArea(BPath, Tol)/(pi));

%% 
% remove inner corners
[BPath_rounded_flipped] = ...
  RemoveAllCorners( FlipPath(BPath), CornerRoundingRadius, Tol, false );
BPath_tmp = FlipPath(BPath_rounded_flipped);

PlotPath(BPath_tmp)

% try different rounding radius before proceeding
BPath = BPath_tmp;

%% 
% remove outer corners
WheelRadius_old = Inf;
WheelRadius_new = (PathPerimeter(BPath,0.00001)/(2*pi))/WheelBezRatio
while abs( WheelRadius_new - WheelRadius_old ) > WheelRadiusTol
  [BPath_new] = ...
    RemoveAllCorners( BPath, WheelRadius_new, Tol, true );
  %
  WheelRadius_old = WheelRadius_new;
  WheelRadius_new = (PathPerimeter(BPath_new,0.00001)/(2*pi))/WheelBezRatio
end
%BPath = BPath_new;
WheelRadius  = WheelRadius_new;
MarkerRadius = WheelRadius*WheelMarkerRatio;
BPath_new = ShiftPath( BPath_new, Shift, Halfen );

if false
  BPath_new = ShiftPath( BPath_new, Shift, Halfen );
  WheelRadius = (PathPerimeter(BPath_new,0.00001)/(2*pi))/WheelBezRatio
  MarkerRadius = WheelRadius*WheelMarkerRatio;
end

%% 
% adjust start point after rounding
PlotPath(BPath_new)

BPath_new = ShiftPath( BPath_new, 2, false);

%%
% colors

ColorVector = {'yellow','magenta', 'red', 'red'};
%ColorVector = {'magenta', 'yellow', 'red', 'red'};

%ColorVector = {'yellow','white', 'red', 'red'};    

%ColorVector = {'yellow','yellow','red', 'red', 'red'};
%ColorVector = {'yellow','red','yellow', 'red', 'red'};

%ColorVector = {'red','white', 'red'};

%%
% preview curve

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;

MarkerAngle0Array = 0;
nPts = size(MarkerAngle0Array,2);

% compute curves
[ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      MaxSpins, CurveOpts);

% plotting per se
figure()
hold on
axis equal
grid on
fill(DecorativeBez(1,:),DecorativeBez(2,:), .15*[1,1,1], 'EdgeColor', 'none')
set(gca,'color', 'k');
for i = 1:nPts
  plot(AllMarkerPos{i}(1,:),AllMarkerPos{i}(2,:),ColorVector{i}, 'LineWidth',2)
end
for i = 1:nPts
  scatter(AllMarkerPos{i}(1,1),AllMarkerPos{i}(2,1),'red','filled','o')
end

%%
% video

% curve parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;

MarkerAngle0Array = 0;
nPts = size(MarkerAngle0Array,2);

% compute curves
[ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      MaxSpins, CurveOpts);

% video parameters
ExtraOpts = {};
ExtraOpts.Plot2Circles = false;
ExtraOpts.Format = 'mp4';
ExtraOpts.Orientation = 'in';
ExtraOpts.Ratio = 16/9;
ExtraOpts.TimerefCurve = 'Average';
ExtraOpts.LineWidth = 2;
ExtraOpts.Tol = Tol;

WhoIsCenter = 1;

% video
MakeVideo_Npts( nPts, WhoIsCenter, WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  40, 10, 'test_251021_20_1', ExtraOpts )
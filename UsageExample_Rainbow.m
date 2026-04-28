% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','LetterC'));
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
BPath_pack = LoadSVG( './curves_svg/octopi1.svg' );
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
BPath = RotatePath( BPath, pi/2 );

% change orientation
BPath = FlipPath(BPath);

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

WheelBezRatio = 12/5;
WheelMarkerRatio = 4/5;

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



% don't remove outer corners
BPath_new = ShiftPath( BPath, Shift, Halfen );
WheelRadius = (PathPerimeter(BPath_new,0.00001)/(2*pi))/WheelBezRatio
MarkerRadius = WheelRadius*WheelMarkerRatio;

%% 
% adjust start point after rounding
PlotPath(BPath_new)

BPath_new = ShiftPath( BPath_new, 1, true);

%%

% declaring common colors
Cmarigold  = [234, 162,  33]/255;
Cwhite     = [255, 255, 255]/255;
Cdkmagenta = [139,   0, 139]/255;

% specific choice
ColorVector = {{Cdkmagenta, Cwhite}, {Cmarigold, Cwhite}};

%%
% preview curve

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 0;

ColorCycles   = 12;
%ColorRefCurve = 'CumDist';
ColorRefCurve = 'Bezier';

aang = 2*pi*(0:1/1:1)+1*pi;
aang(end) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

% compute curves
[ DecorativeBez, AllBezierPos, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

%
color0 = ColorVector{1}{1};
colorF = ColorVector{1}{2};

ColorFunc = cell(1,nPts);
for p = 1:nPts
  switch ColorRefCurve
    case 'CumDist'
      ColorVal = cumsum([0, vecnorm( diff( AllMarkerPos{p}, 1,2), 2,1 )]);
    case 'Bezier'
      ColorVal = cumsum([0, vecnorm( diff( AllBezierPos{p}, 1,2), 2,1 )]);
  end
  ColorVal = ColorVal/ColorVal(end);
  ColorNum = 0.5 + 0.5*cos( ColorCycles* ColorVal * 2*pi );
  ColorFunc{p} = color0' + (colorF'-color0')*ColorNum;
end

% plotting per se
figure();
hold on
axis equal
axis on
grid on
set(gca,'Color','k')

fill(DecorativeBez(1,:),DecorativeBez(2,:), .15*[1,1,1], 'EdgeColor', 'none')

for p = 1:nPts
  drawPts = [ AllMarkerPos{p}, NaN(2,1)];
  drawCol = [ ColorFunc{p}'; NaN(1,3)];
  fill(drawPts(1,:),drawPts(2,:),[0,0,0],'FaceVertexCData',drawCol,'EdgeColor','interp','LineWidth',2)
end

for p = 1:nPts
  switch ColorRefCurve
    case 'CumDist'
      ColorVal = cumsum([0, vecnorm( diff( AllMarkerPos{p}, 1,2), 2,1 )]);
    case 'Bezier'
      ColorVal = cumsum([0, vecnorm( diff( AllBezierPos{p}, 1,2), 2,1 )]);
  end
  ColorVal = ColorVal/ColorVal(end);
  for rref = 0:(0.5/ColorCycles):1
    [~,iidx] = max(ColorVal(ColorVal<=rref));
    scatter(AllMarkerPos{p}(1,iidx),AllMarkerPos{p}(2,iidx),'red','filled','o')
  end
end

%%
% video

% this is a collection of hand-picked colors
%NiceColors = {[255, 59, 209]/255,[165, 36, 61]/255, [208, 241, 191]/255, [240, 45, 58]/255};
%ColorVector = { NiceColors{randi(size(NiceColors,2))} };

% curve parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 1;

MarkerAngle0Array = 0;
nPts = size(MarkerAngle0Array,2);

% compute curves
[ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

% video parameters
ExtraOpts = {};
ExtraOpts.Plot2Circles = false;
ExtraOpts.Format = 'mp4';
ExtraOpts.Orientation = 'in';
ExtraOpts.Ratio = 16/9;
ExtraOpts.TimeRefCurve = 'Average';
%ExtraOpts.TimerefCurve = 'Wheel';
ExtraOpts.LineWidth = 2;
ExtraOpts.Tol = Tol;

ExtraOpts.ColorCycles   = 12;
ExtraOpts.ColorRefCurve = 'CumDist';
%ExtraOpts.ColorRefCurve = 'Bezier';

WhoIsCenter = 1;

% video
MakeVideo_Npts( nPts, WhoIsCenter, WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  30, 7.5, 'test_260428_21_1', ExtraOpts )
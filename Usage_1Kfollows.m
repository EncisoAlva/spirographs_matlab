% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','OG_example'));
BPath = BPath_pack{1};
clear BPath_pack

%%
% show control points

BPath = ShiftPath( BPath, -1, true );
BPath = FlipPath(BPath);

PlotPath(BPath)

%%
% parameters

% technical stuff
Tol = 0.005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.000001;

% designer stuff
MarkerAngle0 = 0;

%WheelBezRatio = 5/3;
WheelBezRatio = 24/5;
WheelMarkerRatio = 4/5;

Shift  = 0;
Halfen = false;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.001*PathArea(BPath, Tol)/(pi));

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

%%
% colors

PlotPath(BPath_new)

ColorVector = {'yellow'};

%%
% preview curve

% MarkerRadius = WheelRadius;

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 0;

aang = 2*pi*(0:1/2:1)+0*pi;
aang(end) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

% compute curves
[ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

% plotting per se
figure()
hold on
axis equal
grid on
fill(DecorativeBez(1,:),DecorativeBez(2,:), .15*[1,1,1], 'EdgeColor', 'none')
set(gca,'color', 'k');
for i = 1:nPts
  plot(AllMarkerPos{i}(1,:),AllMarkerPos{i}(2,:),'Color',ColorVector{mod(i-1,k)+1}, 'LineWidth',2)
end
for i = 1:nPts
  scatter(AllMarkerPos{i}(1,1),AllMarkerPos{i}(2,1),'red','filled','o')
end

%%

[ ~,...
  AllBezierPos1, AllLocTime1, ...
  AllWhCtrPos1, AllMarkerPos1, AllMarkerAngle1 ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

WheelRadii1 = WheelRadius*ones(1,nPts);

%%
% parameters

WheelBezRatio = 24/5;

% don't remove outer corners
BPath_new = ShiftPath( BPath_new, Shift, Halfen );
Perimeter = (PathPerimeter(BPath_new,0.00001)/(2*pi));
WheelRadius  = Perimeter/WheelBezRatio

%%
% colors

BPath_new = FlipPath(BPath_new);

ColorVector = {'magenta'};

%%
% preview curve

% MarkerRadius = WheelRadius;

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 0;

aang = 2*pi*(0:1/2:1)+1*pi;
aang(end) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

% compute curves
[ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

% plotting per se
figure()
hold on
axis equal
grid on
fill(DecorativeBez(1,:),DecorativeBez(2,:), .15*[1,1,1], 'EdgeColor', 'none')
set(gca,'color', 'k');
for i = 1:nPts
  plot(AllMarkerPos{i}(1,:),AllMarkerPos{i}(2,:),'Color',ColorVector{mod(i-1,k)+1}, 'LineWidth',2)
end
for i = 1:nPts
  scatter(AllMarkerPos{i}(1,1),AllMarkerPos{i}(2,1),'red','filled','o')
end

%%

CurveOpts.ChangeOrient = [true, true, true];

[ DecorativeBez,...
  AllBezierPos2, AllLocTime2, ...
  AllWhCtrPos2, AllMarkerPos2, AllMarkerAngle2 ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

WheelRadii2 = WheelRadius*ones(1,nPts);

%%

ColorVector = {'yellow', 'yellow', 'magenta', 'magenta'};

AllBezierPos   = [AllBezierPos1,   AllBezierPos2];
AllLocTime     = [AllLocTime1,     AllLocTime2];
AllWhCtrPos    = [AllWhCtrPos1,    AllWhCtrPos2];
AllMarkerPos   = [AllMarkerPos1,   AllMarkerPos2];
AllMarkerAngle = [AllMarkerAngle1, AllMarkerAngle2];

%%
% video

% video parameters
ExtraOpts = {};
ExtraOpts.Plot2Circles = false;
ExtraOpts.Format = 'mp4';
ExtraOpts.Orientation = 'in';
ExtraOpts.Ratio = 16/9;
%ExtraOpts.TimerefCurve = 'Average';
%ExtraOpts.TimerefCurve = 'Wheel';
ExtraOpts.TimeRefCurve = 'Bezier';
ExtraOpts.LineWidth = 2;
ExtraOpts.Tol = Tol;
ExtraOpts.WheelRadii = [WheelRadii1, WheelRadii2];
%
ExtraOpts.WhoIsCenter = [1,4];

% video
MakeVideo_Npts( 4, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  30, 7.5, 'test_1K_3', ExtraOpts )
% Rolling wheel has a convex hole, and he marker slides throught the hole.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load path to roll on

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','TwoCircles'));
BPath = BPath_pack{1};

clear BPath_pack

%%
% check available curves in the example file
who -file ExampleCollections.mat

% load curve
BPath_pack1 = struct2cell(load('ExampleCollections.mat','Circlegon'));
BPath_pack2 = BPath_pack1{1};
BPath = BPath_pack2{2};

clear BPath_pack1 BPath_pack2

%%
% load from file
BPath_pack = LoadSVG( './curves_svg/Yscavenge.svg' );
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

BPath = ShiftPath( BPath, 1, false );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load path for the hole in the wheel

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','Heart'));
HPath = BPath_pack{1};

clear BPath_pack

%%
% check available curves in the example file
who -file ExampleCollections.mat

% load curve
BPath_pack1 = struct2cell(load('ExampleCollections.mat','Circlegon'));
BPath_pack2 = BPath_pack1{1};
HPath = BPath_pack2{2};

clear BPath_pack1 BPath_pack2

%%
% load from file
BPath_pack = LoadSVG( './curves_svg/Yscavenge.svg' );
HPath = BPath_pack{1};

clear BPath_pack

%%
% pre-processing
HPath = RemovePointCurves( HPath, 0.0001 );

% this one is because I'm using absolute tolerance instead of relative
HPath = RescalePath( HPath, 0.5, 0.5 );

% line with bad encoding, the normal vector will be wrong
HPath = ForceCubicLines( HPath );

% rotate half a spin
for i = 1:size(HPath, 2)
  HPath{i} = [1,0; 0,-1] * HPath{i};
end

% rotate by an angle
HPath = RotatePath( HPath, -pi/2 );

% move by a vector
HPath = TranslatePath( HPath, [0,0.75]' );

% change orientation
HPath = FlipPath(HPath);

%%
% show control points

PlotHole(HPath)

HPath = ShiftPath( HPath, -1, false );

%%
% prepare for interpolation
SetupHole(HPath, 0.001, true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% parameters

% technical stuff
Tol = 0.005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.000001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 6/5;

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
%MarkerRadius = WheelRadius*WheelMarkerRatio;
BPath_new = ShiftPath( BPath_new, Shift, Halfen );



% don't remove outer corners
BPath_new = ShiftPath( BPath, Shift, Halfen );
WheelRadius = (PathPerimeter(BPath_new,0.00001)/(2*pi))/WheelBezRatio
%MarkerRadius = WheelRadius*WheelMarkerRatio;

%% 
% adjust start point after rounding
PlotPath(BPath_new)

BPath_new = ShiftPath( BPath_new, 1, true);

%%
% colors

%ColorVector = {'yellow','magenta','blue','red','green'};

%ColorVector = {'red','white', 'yellow'};

ColorVector = {'white', 'red'};

%ColorVector = {'white'};

%ColorVector = {'red'};

%ColorVector = {'yellow'};

%ColorVector = {'magenta'};

%ColorVector = {'green'};

%ColorVector = {'magenta','yellow'};

%ColorVector = {'yellow', 'magenta', 'magenta'};

%%
% preview curve

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 1;

aang = 2*pi*(0:1/1:1)+pi*1;
aang(end) = [];
%aang(1) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

% compute curves
[ DecorativeBez, AllBezierPos, ~, AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
    SetupCurves_Npts( 1, BPath_new, WheelRadius, WheelRadius, MarkerAngle0Array, ...
      CurveOpts);

% now computing the marker sliding through the hole
for i = 1:nPts
  cumDist  = cumsum([0,vecnorm(diff(AllBezierPos{i},1,2),2,1)]);
  locAngle = mod( cumDist * PathPerimeter(BPath,Tol) / ( WheelRadius * 2*pi ), 2*pi );
  %
  MarkerPos_tmp = zeros(size(AllWhCtrPos{i}));
  MarkerPos_tmp(1,:) = WheelRadius*interp1(TH_vector, H_vector(1,:), mod(locAngle,2*pi));
  MarkerPos_tmp(2,:) = WheelRadius*interp1(TH_vector, H_vector(2,:), mod(locAngle,2*pi));
  %
  for j = 1:size(AllWhCtrPos{i},2)
    th = AllMarkerAngle{i}(j) ;
    %
    AllMarkerPos{i}(:,j) = AllWhCtrPos{i}(:,j) +[cos(th) -sin(th); sin(th) cos(th)]* MarkerPos_tmp(:,j);
  end
end

figure()
plot(locAngle)
figure()
plot(AllMarkerAngle{i})

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
% fancy plot, for ads

%ColorVector = {'white'};

%ColorVector = {'red'};

%ColorVector = {'yellow'};

%ColorVector = {'magenta'};

%ColorVector = {'green'};

%ColorVector = {'cyan'};

f1 = figure('Name','Just the curve');
hold on
axis equal
axis off
set(gcf,'Color','k')
set(gca,'Color','k')

plot(DecorativeBez(1,:),DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)

for p = 1:nPts
  plot(AllMarkerPos{p}(1,:),AllMarkerPos{p}(2,:),'Color',ColorVector{p},'LineWidth',2)
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
CurveOpts.MaxSpins = 1;
CurveOpts.MinSpins = 0;

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
ExtraOpts.TimerefCurve = 'Average';
%ExtraOpts.TimerefCurve = 'Wheel';
ExtraOpts.LineWidth = 2;
ExtraOpts.Tol = Tol;

WhoIsCenter = 1;

% video
MakeVideo_Npts( nPts, WhoIsCenter, WheelRadius, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  30, 7.5, 'test_260324_03_3', ExtraOpts )
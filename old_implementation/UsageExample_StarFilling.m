% One circle rolling on the outside, marking 2 points. Optional is to round
% the corners.
%

%%
% check available curves in the example file
who -file ExampleCurves.mat

% load curve
BPath_pack = struct2cell(load('ExampleCurves.mat','Guinivere'));
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
Halfen = true;

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
% colors

%ColorVector = {'white'};

%ColorVector = {'red'};

%ColorVector = {'yellow'};

%ColorVector = {'magenta'};

%ColorVector = {'magenta','yellow'};

%ColorVector = {'yellow', 'magenta', 'magenta'};


%%
% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 0;

aang = 2*pi*(0:1/1:1);
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
% generate points according to the following algorithm:
% 1. create points uniformly on a square that coers the curve
% 2. consider only points inside the curve
% 3. compute the distance curve-point, if the point is far from the curve,
%    then it is more likely to be removed
% 4. comb the points based on its distance to the curve

% parameters
nStars = 5000;
starRadius = 0.001;

% create box containing the curve
Xmin =  Inf;
Xmax = -Inf;
Ymin =  Inf;
Ymax = -Inf;
for i = 1:nPts
  Xmin = min( Xmin, min(AllMarkerPos{i}(1,:)) );
  Xmax = max( Xmax, max(AllMarkerPos{i}(1,:)) );
  Ymin = min( Ymin, min(AllMarkerPos{i}(2,:)) );
  Ymax = max( Ymax, max(AllMarkerPos{i}(2,:)) );
end

% main loop
rng(2026)
WB = waitbar(0,strcat('Generating points...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2026)');
%
coordStar = zeros(2, nStars);
currStar = 1;
while currStar < nStars
  %
  % 
  if getappdata(WB,'canceling')
    disp('Ended by user.')
    close(v);
    delete(WB)
    break
  end
  %
  xp = Xmin + (Xmax - Xmin) * rand(1, 1);
  yp = Ymin + (Ymax - Ymin) * rand(1, 1);
  %
  % if point in curve
  if inpolygon(xp, yp, AllMarkerPos{1}(1,:), AllMarkerPos{1}(2,:))
    distStar = min( vecnorm( AllMarkerPos{1} - [xp; yp], 2, 1 ) );
    probKeep = 1 - distStar/(distStar+starRadius);
    %
    % trimming
    if rand(1, 1) <= probKeep
      coordStar(:,currStar) = [xp; yp];
      %
      % update progressbar
      waitbar(currStar/nStars,WB);
      currStar = currStar + 1;
    end
  end
end

% delete waitbar, if it is still around
if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% check
figure()
hold on
axis equal
grid on
set(gca,'color', 'k');
for i = 1:nPts
  plot(AllMarkerPos{i}(1,:),AllMarkerPos{i}(2,:),'Color',ColorVector{mod(i-1,k)+1}, 'LineWidth',2)
end
scatter(coordStar(1,:),coordStar(2,:),3,[48, 133, 194]/255,'filled','o')

% fancy options
x0 = min(DecorativeBez(1,:));
xF = max(DecorativeBez(1,:));
y0 = min(DecorativeBez(2,:));
yF = max(DecorativeBez(2,:));
for p = 1:nPts
  x0 = min( x0, min(AllMarkerPos{p}(1,:)) );
  xF = max( xF, max(AllMarkerPos{p}(1,:)) );
  y0 = min( y0, min(AllMarkerPos{p}(2,:)) );
  yF = max( yF, max(AllMarkerPos{p}(2,:)) );
end
%
x_ran = xF - x0;
y_ran = yF - y0;
%
if y_ran/x_ran < ExpectedRatio
  y_ran_new = x_ran*ExpectedRatio;
  y0 = y0 - (y_ran_new - y_ran)/2;
  yF = yF + (y_ran_new - y_ran)/2;
else 
  if y_ran/x_ran > ExpectedRatio
    x_ran_new = y_ran/ExpectedRatio;
    x0 = x0 - (x_ran_new - x_ran)/2;
    xF = xF + (x_ran_new - x_ran)/2;
  end
end
%
xlim([x0 xF])
ylim([y0 yF])

% the actual figure
figure()
hold on
axis equal
axis off
grid off
set(gca,'color', 'k');
scatter(coordStar(1,:),coordStar(2,:),2,[48, 133, 194]/255,'filled','o')


xlim([x0 xF])
ylim([y0 yF])

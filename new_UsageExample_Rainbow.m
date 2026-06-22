%%
% check available curves
BezPath.CheckExamples();

%%
% generate empty curve
Curve = BezGlissette();

% load path from example
Curve.LoadBezPath( 'UniqueCurve', 'LetterC' )

% load path from indexed example
Curve.LoadBezPath( 'IndexedCurve', 'Circlegon', 2 )

% load path from SVG file
Curve.LoadBezPath( 'SVG', './curves_svg/Yscavenge.svg' )

%%
% pre-processing

Curve.BPath.Rotate( pi/2 );

Curve.BPath.Flip();

Curve.BPath.Shift( 1, false );

Curve.BPath.PlotPath()

%%
% designer parameters

% designer stuff
WheelBezRatio = 6;
WheelMarkerRatio = 4/5;

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

% declaring common colors
Cmarigold  = [234, 162,  33]/255;
Cwhite     = [255, 255, 255]/255;
Cdkmagenta = [139,   0, 139]/255;
Cscarlet   = [255,  36,   0]/255;
Cblack     = [  0,   0,   0];
Cruby_cl   = [224,  17,  95]/255;
Cruby_dk   = [ 78,   9,  15]/255;
Cbatman_y  = [152, 136,  41]/255;
Cpurple    = [128,   0, 128]/255;
Cpurple_cl = [192,   0, 192]/255;
Cblue_fcy  = [ 36, 122, 253]/255;
Cflower_y  = [255, 229,  90]/255;

% specific choice
%ColorVector = {{Cwhite, Cmarigold, Cmarigold}, {Cscarlet, Cwhite}};

ColorVector = {{Cblue_fcy, Cblue_fcy, Cflower_y, Cpurple_cl, Cpurple_cl, Cflower_y}};

%%
% preview curve

% setup parameters
CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;
CurveOpts.MaxSpins = 100;
CurveOpts.MinSpins = 6;

ColorCycles   = 3;
ColorRefCurve = 'CumDist';
%ColorRefCurve = 'Bezier';

aang = 2*pi*(0:1/2:1)+0*pi;
aang(end) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

% compute curves
[ DecorativeBez, AllBezierPos, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

ColorFunc = cell(1,nPts);
for p = 1:nPts
  switch ColorRefCurve
    case 'CumDist'
      ColorVal = cumsum([0, vecnorm( diff( AllMarkerPos{p}, 1,2), 2,1 )]);
    case 'Bezier'
      ColorVal = cumsum([0, vecnorm( diff( AllBezierPos{p}, 1,2), 2,1 )]);
  end
  ColorVal = ColorVal/ColorVal(end);
  ColorNum = 0.5 - 0.5*cos( ColorCycles* ColorVal * 2*pi );
  %
  colorTable = zeros(size(ColorVector{p},2), 3+1);
  for q = 1:size(ColorVector{p},2)
    colorTable(q,2:end) = ColorVector{p}{q};
  end
  colorTable(:,1) = linspace(0,1, size(ColorVector{p},2));
  %
  ColorFunc{p} = zeros(3, size(ColorNum,2));
  ColorFunc{p}(1,:) = interp1(colorTable(:,1),colorTable(:,2), ColorNum);
  ColorFunc{p}(2,:) = interp1(colorTable(:,1),colorTable(:,3), ColorNum);
  ColorFunc{p}(3,:) = interp1(colorTable(:,1),colorTable(:,4), ColorNum);
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
ExtraOpts.WheelRadii = WheelRadius*ones(1,nPts);

ExtraOpts.ColorCycles   = 3;
ExtraOpts.ColorRefCurve = 'CumDist';
%ExtraOpts.ColorRefCurve = 'Bezier';
%
ExtraOpts.WhoIsCenter = 1;

% video
MakeVideo_Npts( nPts, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  30, 7.5, 'test_2600515_02_3', ExtraOpts )
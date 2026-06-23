%%
% check available curves
BezPath.CheckExamples();

%%
% generate empty curve
Curve = BezGlissette( 'Standard' );

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
Curve.Set_Wheel1BezRatio( 6 );

Curve.WheelMarkerRatio = 4/5;

%%
% additional design parameters
Curve.RemoveCorners_Rolling = true;
Curve.RemoveCorners_NonRolling = true;

%% 
% pre-processing (optional)

Curve.RemoveCorners()
Curve.BPath.PlotPath()

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

Curve.SetColor( {Cblue_fcy, Cblue_fcy, Cflower_y, Cpurple_cl, Cpurple_cl, Cflower_y},...
 'CumDist', 3 );
%'CumDist'
%'Bezier'

%%
% preview curve

% compute curves
[ DecorativeBez, AllBezierPos, ~, ~, AllMarkerPos, ~ ] = ...
    SetupCurves_Npts( nPts, BPath_new, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
      CurveOpts);

Curve.PlotGlissette()

%%
% video

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
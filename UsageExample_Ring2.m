% The circular sliding wheel has a circular hole inside of it, and thus it
% is referred to as a ring. Inside the ring, there is a smaller circular
% wheel rolling; the marker is in this second wheel.

%%
% check available curves
BezPath.CheckExamples();

%%
% generate empty curve
Curve = BezGlissette( 'Ring2' );

% load path from example
%Curve.LoadRollingPath( 'UniqueCurve', 'LetterC' )

% load path from indexed example
Curve.LoadRollingPath( 'IndexedCurve', 'Circlegon', 2 )

% load path from SVG file
%Curve.LoadRollingPath( 'SVG', './curves_svg/Yscavenge.svg' )

%%
% pre-processing

%Curve.BPath.Rotate( pi/2 );

Curve.BPath.Flip();

Curve.BPath.Shift( 2 );

Curve.PlotPath()

%%
% design parameters
%Curve.RemoveCorners_Rolling = true;
%Curve.RemoveCorners_NonRolling = true;

% update all relevant parameter with each change of ratios
Curve.AutoUpdate = true;

%%
% designer parameters

%>>>
WheelBezRatio  = 120/96;
HoleBezRatio   = 120/40;
Wheel2BezRatio = 120/32;
Wheel2Marker2Ratio = 4/5;
%<<<

Curve.Set_Wheel1BezRatio( 6 );

%>>>
Perimeter = (PathPerimeter(BPath_new,0.00001)/(2*pi));
WheelRadius  = WheelRadius_new;
Wheel2Radius = Perimeter/Wheel2BezRatio;
HoleRadius   = Perimeter/HoleBezRatio;
Marker2Radius = Wheel2Radius*Wheel2Marker2Ratio;
BPath_new = ShiftPath( BPath_new, Shift, Halfen );
%<<<

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

Curve.SetColor( {Cflower_y, Cblue_fcy},...
 'CumDist', 6 );
%'CumDist'
%'Bezier'

%%
% preview curve

Curve.ProcessGlissette()

Curve.PlotGlissette()

%%
% fast change

Curve.Set_Wheel1BezRatio( 5,3 );

Curve.ProcessGlissette()

%Curve.ColorCycles = 5;
Curve.PlotGlissette()

%%
% video

% video parameters
VideoOpts = {};
VideoOpts.Format = 'mp4';
VideoOpts.Orientation = 'in';
VideoOpts.Ratio = 16/9;
VideoOpts.LineWidth = 2;
VideoOpts.WhoIsCenter = 1;
VideoOpts.WheelRadii = Curve.Wheel1Radius;
%VideoOpts.TimeRefCurve = 'Average';
%VideoOpts.TimeRefCurve = 'Wheel';
%VideoOpts.TimeRefCurve = 'Marker';
VideoOpts.TimeRefCurve = 'Avg_MarkerBezier';

VideoOpts.AddDateTimeIndex = true;

% video
Curve.MakeVideo( 'test', VideoOpts )

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
CurveOpts.MinSpins = 3;

MarkerAngle0Array = 0;
nPts = size(MarkerAngle0Array,2);

% compute curves
[ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_Ring2( BPath_new, WheelRadius, ...
    Wheel2Radius, Marker2Radius, HoleRadius, CtrHoleDist, ...
    MarkerAngle0Array, ...
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
ExtraOpts.WheelRadii = WheelRadius*ones(1,nPts);
%
ExtraOpts.WhoIsCenter = 1;

% video
MakeVideo_Npts( nPts, ...
  DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle,...
  ColorVector, ...
  30, 7.5, 'test_260409_14_1', ExtraOpts )
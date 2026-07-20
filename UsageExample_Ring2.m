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

% the parameters for the ring are more complex
% for compatibility with the terminology of Wild Gears, the sizes of the
% circles are given in terms of their perimeters as integers
WheelPerimeters = [];
WheelPerimeters.RollingCurve = 96;
WheelPerimeters.Wheel1_Outer = 80;
WheelPerimeters.Wheel1_Inner = 40;
WheelPerimeters.Wheel2       = 32;
WheelPerimeters.Wheel2MarkerRatio = 1;

Curve.Set_RingParameters( WheelPerimeters );

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
 'CumDist', 1 );
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
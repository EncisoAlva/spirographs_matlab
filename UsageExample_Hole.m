%%
% check available curves
BezPath.CheckExamples();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load path to roll on

%%
% generate empty curve
Curve = BezGlissette( 'Hole' );

% load path from example
Curve.LoadRollingPath( 'UniqueCurve', 'LetterC' )

% load path from indexed example
Curve.LoadRollingPath( 'IndexedCurve', 'Circlegon', 2 )

% load path from SVG file
Curve.LoadRollingPath( 'SVG', './curves_svg/Yscavenge.svg' )

%%
% pre-processing

Curve.BPath.Rotate( pi/2 );

Curve.BPath.Flip();

Curve.BPath.Shift( 1, false );

Curve.BPath.PlotPath()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load path for the hole in the wheel

% load path from example
Curve.LoadHolePath( 'UniqueCurve', 'LetterC' )

% load path from indexed example
Curve.LoadHolePath( 'IndexedCurve', 'Circlegon', 2 )

% load path from SVG file
Curve.LoadHolePath( 'SVG', './curves_svg/Yscavenge.svg' )

%%
% pre-processing

Curve.BPath.FitBox( [0;0], [1;1] )

Curve.BPath.Flip().RemoveAllCorners().Flip()

Curve.BPath.Translate( [0,0.5]' )

Curve.BPath.Rotate( -pi/2 );

Curve.BPath.Flip();

Curve.BPath.Shift( 1, false );

Curve.BPath.PlotHole()

%%
% prepare for interpolation
Curve.SetupHole( true )

% experimental
Curve.DEV_SetupHole_concave( true )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% designer parameters
Curve.Set_Wheel1BezRatio( 8, 7 );

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

Curve.ProcessGlissette()

Curve.PlotGlissette()

%%
% video

% video parameters
VideoOpts = {};
VideoOpts.Format = 'mp4';
VideoOpts.Orientation = 'in';
VideoOpts.Ratio = 16/9;
%ExtraOpts.TimeRefCurve = 'Average';
%ExtraOpts.TimeRefCurve = 'Wheel';
%ExtraOpts.TimeRefCurve = 'Marker';
ExtraOpts.TimeRefCurve = 'Avg_MarkerBezier';

VideoOpts.LineWidth = 2;
VideoOpts.WhoIsCenter = 1;
VideoOpts.WheelRadii = Curve.Wheel1Radius;

% video
Curve.MakeVideo( 'test_2600515_02_3', VideoOpts )
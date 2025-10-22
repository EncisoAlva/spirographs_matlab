% Given a list of Bezier curves and parameters, creates arrays 
% necessary to create the animations.
% The goal is to test multiple parameter configurations before making the 
% animation, since testing parameters is much faster than the rest of the
% animation.
%
% ---- INUPUT ------------------------------------------------------------
%    CtrlPtsArray  Array with control points for each one of the Bezier
%                  curves that make the curve {?} <- [2,4]'s
%   WheelBezRatio  Perimeter of shape / Perimeter of wheel [1]
% WheelMarkerRatio Center of wheel to marker / Radius of wheel [1]
%    MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                  the wheelcenter-marker line [1]
%    MaxDistDelta  Maxi allowable distance between neighboring points [1]
%        CloseTol  Max distance between first and last point [1]
%        MaxSpins  Max full rotations of wheel around whole shape [1]
%
% ---- OUTPUT ------------------------------------------------------------
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%   WhCtrPos1/2  Location of the wheel center at timestamps [2x?]
%  MarkerPos1/2  Location of marker at timepoints [2x?]
%MarkerAngle1/2  Angle of marker, at each timepoint, with respect to the
%                x-axis [1x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_4pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    ScaleFactor, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, ExtraOpts)

% this evaluation is for background decoration only
DecorativeBez = AllBezierEval(CtrlPtsArray, MaxDistDelta);

% containers
AllBezierPos   = {};
AllLocTime     = {};
AllWhCtrPos    = {};
AllMarkerPos   = {};
AllMarkerAngle = {};

% outer curves
[ DecorativeBez,...
  AllBezierPosA, AllLocTimeA, ...
  AllWhCtrPosA, AllMarkerPosA, AllMarkerAngleA ] = ...
  SetupCurves_2pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, ExtraOpts);

% inner curves
[ ~,...
  AllBezierPosB, AllLocTimeB, ...
  AllWhCtrPosB, AllMarkerPosB, AllMarkerAngleB ] = ...
  SetupCurves_2pts( FlipBezierAll(CtrlPtsArray), WheelRadius/ScaleFactor, MarkerRadius/ScaleFactor, -MarkerAngle0*ScaleFactor+pi, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins, ExtraOpts);

AllBezierPosB{1}   = flip(AllBezierPosB{1}, 2);
AllWhCtrPosB{1}    = flip( AllWhCtrPosB{1}, 2 );
AllMarkerPosB{1}   = flip(AllMarkerPosB{1}, 2);
AllMarkerAngleB{1} = flip(AllMarkerAngleB{1}, 2);
%
AllBezierPosB{2}   = flip(AllBezierPosB{2}, 2);
AllWhCtrPosB{2}    = flip(AllWhCtrPosB{2}, 2);
AllMarkerPosB{2}   = flip(AllMarkerPosB{2}, 2);
AllMarkerAngleB{2} = flip(AllMarkerAngleB{2}, 2);

AllLocTimeB{1} = max(AllLocTimeB{1}) - flip(AllLocTimeB{1});
AllLocTimeB{2} = max(AllLocTimeB{2}) - flip(AllLocTimeB{2});

% multiplicity due to scaling factor
if ScaleFactor > 1
  AllBezierPosB{1}   = kron(ones(1,ScaleFactor), AllBezierPosB{1});
  AllBezierPosB{2}   = kron(ones(1,ScaleFactor), AllBezierPosB{2});
  AllLocTimeB{1}     = kron(ones(1,ScaleFactor), AllLocTimeB{1}) + kron( AllLocTimeB{1}(end)*(0:(ScaleFactor-1)), ones(size(AllLocTimeB{1})) );
  AllLocTimeB{2}     = kron(ones(1,ScaleFactor), AllLocTimeB{2}) + kron( AllLocTimeB{2}(end)*(0:(ScaleFactor-1)), ones(size(AllLocTimeB{2})) );
  AllWhCtrPosB{1}    = kron(ones(1,ScaleFactor), AllWhCtrPosB{1});
  AllWhCtrPosB{2}    = kron(ones(1,ScaleFactor), AllWhCtrPosB{2});
  AllMarkerPosB{1}   = kron(ones(1,ScaleFactor), AllMarkerPosB{1});
  AllMarkerPosB{2}   = kron(ones(1,ScaleFactor), AllMarkerPosB{2});
  AllMarkerAngleB{1} = kron(ones(1,ScaleFactor), AllMarkerAngleB{1});
  AllMarkerAngleB{2} = kron(ones(1,ScaleFactor), AllMarkerAngleB{2});
end

% Store results in containers
AllBezierPos = [ AllBezierPosA, AllBezierPosB ];
AllLocTime   = [ AllLocTimeA, AllLocTimeB ];
AllWhCtrPos  = [ AllWhCtrPosA, AllWhCtrPosB ];
AllMarkerPos = [ AllMarkerPosA, AllMarkerPosB ];
AllMarkerAngle = [ AllMarkerAngleA, AllMarkerAngleB ];

end
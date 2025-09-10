% Given a list of Bezier curves and parameters, creates arrays 
% necessary to create the animations.
% The goal is to test multiple parameter configurations before making the 
% animation, since testing parameters is much faster than the rest of the
% animation.
%
% ---- INUPUT ------------------------------------------------------------
%  CtrlPtsArray  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%  MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                wheelcenter-marker line [1]
%  MaxDistDelta  Maximum allowable distance between neighboring points [1]
%      CloseTol  Max distance between first and last point [1]
%      MaxSpins  Max full rotations of wheel around whole shape [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      WhCtrPos  Location of the wheel center at timestamps [2x?]
%     MarkerPos  Location of marker at timepoints [2x?]
%   MarkerAngle  Angle of marker, at each timepoint, with respect to the
%                x-axis [1x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [Time, WhCtrPos, MarkerPos, MarkerAngle] = ...
  SetupCurve_2pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins)

% specific to this example
WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/(15/7);
MarkerRadius = WheelRadius*(1);

% Bezier curve
BezierPos = AllBezierEval( CtrlPtsArray, MaxDistDelta );

% spirograph curve
[~, WhCtrPos1, MarkerPos1, MarkerAngle1] = ...
  AllBeziers( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

[~, WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  AllBeziers( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0+pi, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins);

% time, parametrized by the arc of the marker or the ar of the wheel center
TimeFromMarker = zeros(1,size(MarkerPos1,2));
TimeFromMarker(2:end) = cumsum( vecnorm( diff(MarkerPos1,1,2), 2, 1 ) );

TimeFromWheel1 = zeros(1,size(WhCtrPos1,2));
TimeFromWheel1(2:end) = cumsum( vecnorm( diff(WhCtrPos1,1,2), 2, 1 ) );

TimeFromWheel2 = zeros(1,size(WhCtrPos2,2));
TimeFromWheel2(2:end) = cumsum( vecnorm( diff(WhCtrPos2,1,2), 2, 1 ) );

% patch
MarkerPos1(:,end+1) = MarkerPos1(:,1);
MarkerAngle1(end+1) = MarkerAngle1(1);

MarkerPos2(:,end+1) = MarkerPos2(:,1);
MarkerAngle2(end+1) = MarkerAngle2(1);

end
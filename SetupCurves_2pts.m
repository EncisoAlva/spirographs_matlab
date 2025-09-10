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
function [WheelRadius, MarkerRadius, BezierPos, ...
  WhCtrPos1, MarkerPos1, MarkerAngle1,...
  WhCtrPos2, MarkerPos2, MarkerAngle2] = ...
  SetupCurves_2pts( CtrlPtsArray, WheelBezRatio, WheelMarkerRatio, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins)

% specific to this example
WheelRadius  = (BezierPerimeter(CtrlPtsArray,0.00001)/(2*pi))/WheelBezRatio;
MarkerRadius = WheelRadius*WheelMarkerRatio;

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

% patch
MarkerPos1(:,end+1) = MarkerPos1(:,1);
MarkerAngle1(end+1) = MarkerAngle1(1);

MarkerPos2(:,end+1) = MarkerPos2(:,1);
MarkerAngle2(end+1) = MarkerAngle2(1);

end
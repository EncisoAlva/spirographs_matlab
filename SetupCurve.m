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
  SetupCurve( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins)
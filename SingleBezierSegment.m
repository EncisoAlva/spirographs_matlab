% Create the spirograph for one Bezier cubic curve; the full curve is the
% concatenation of multiple Bezier curves.
%
% ---- INUPUT ------------------------------------------------------------
%       CtrlPts  Control points for Bezier curve, [2x4]
%     CirRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%  MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                wheelcenter-marker line [1]
%    MarkerPos0  Initial position of the marker [2x1]
%       CirPos0  Initial position of the wheel center [2x1]
%         Time0  Initial timestamp [1]
%  MaxDistDelta  Maximum allowable distance between neighboring points [1]
%   MarkerSpeed  Speed at which the marker goes throuh its arc [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      CirCentT  Location of the wheel center at timestamps [2x?]
%       MarkerT  Location of marker at timepoints [2x?]
%
% *Notice that some parameters are redundant; this is to avoid computing
% multile times the same parameters.
%
function [Time, CirCenterT, MarkerT] = ...
  SingleBezierSegment( CtrlPts, CirRadius, MarkerRadius, ...
  MarkerAngle0, MarkerPos0, CirPos0, Time0, MaxDistDelta, MarkerSpeed )

LocalTime = 0:(1/MaxDistDelta):1;

end
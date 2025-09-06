% Special treatment for the 'corner' case: there is a non-zero angle
% between two contiguous Bezier curves. In a real spirograph, this means
% that the wheel will rotate without rolling, keeping the contact point
% with the curve but changing its center.
% On a more practical sense, the center of the wheel will follow a circle
% arc, whose center is the contact point with the Bezier curves (which is a
% common point for both curves) and whose radius is the wheel radius. The
% starting and ending points are given by the tangents to the Bezier
% curves.
%
% ---- INUPUT ------------------------------------------------------------
%       CtrlPts  Control points for Bezier curve, [2x4]
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%  MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                wheelcenter-marker line [1]
%         Time0  Initial timestamp [1]
%  MaxDistDelta  Maximum allowable distance between neighboring points [1]
%   MarkerSpeed  Speed at which the marker goes throuh its arc [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      WhCtrPos  Location of the wheel center at timestamps [2x?]
%     MarkerPos  Location of marker at timepoints [2x?]
%
% *Notice that some parameters are redundant; this is to avoid computing
% multile times the same parameters.
%
function [Time, WhCtrPosi, MarkerPos, MarkerAngle] = ...
  RollCorner( CtrlPtsPre, CtrlPtsPos, WheelRadius, MarkerRadius, ...
  MarkerAngle0, Time0, MaxDistDelta, MarkerSpeed )

% compute arc angle that the wheel center will describe
WhNormalPre = EvalBezierNormal( CtrlPtsPre, 1, WheelRadius );
WhNormalPos = EvalBezierNormal( CtrlPtsPos, 0, WheelRadius );
WhCtrPre = CtrlPtsPre(:,end) + WhNormalPre;
CornerAngle = atan2(WhNormalPos(2),WhNormalPos(1)) - atan2(WhNormalPre(2),WhNormalPre(1));

% compute the length that the marker will describe
MarkerPos0 = WhCtrPre + [cos(MarkerAngle0); sin(MarkerAngle0)]*MarkerRadius;
MarkerArcLength = abs(CornerAngle) * norm( MarkerPos0 - CtrlPtsPos(:,1) );

% adjust the rotation, time is set in terms of the marker
LocalTime   = 0:1/ceil(MarkerArcLength / MaxDistDelta):1;
MarkerAngle = MarkerAngle0 + LocalTime*CornerAngle;

LocalWhCtrAngle   = atan2(WhNormalPre(2),WhNormalPre(1));
LocalMarkerAngle  = atan2(MarkerPos0(2)-CtrlPtsPre(2,end),MarkerPos0(1)-CtrlPtsPre(1,end));
LocalMarkerRadius = norm( MarkerPos0 - CtrlPtsPos(:,1) );

WhCtrPosi = CtrlPtsPre(:,end) + [cos(LocalWhCtrAngle  + LocalTime*CornerAngle); sin(LocalWhCtrAngle  + LocalTime*CornerAngle)]*WheelRadius;
MarkerPos = CtrlPtsPre(:,end) + [cos(LocalMarkerAngle + LocalTime*CornerAngle); sin(LocalMarkerAngle + LocalTime*CornerAngle)]*LocalMarkerRadius;

Time = LocalTime*MarkerSpeed + Time0;

% technical
MarkerAngle = mod(MarkerAngle, 2*pi);

end

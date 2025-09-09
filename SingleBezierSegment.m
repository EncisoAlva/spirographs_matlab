% Create the spirograph for one Bezier cubic curve; the full curve is the
% concatenation of multiple Bezier curves.
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
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      WhCtrPos  Location of the wheel center at timestamps [2x?]
%     MarkerPos  Location of marker at timepoints [2x?]
%   MarkerAngle  Angle of marker, at each timepoint, with respect to the
%                x-axis [1x?]
%
% *Notice that some parameters are redundant; this is to avoid computing
% multile times the same parameters.
%
function [Time, WhCtrPos, MarkerPos, MarkerAngle] = ...
  SingleBezierSegment( CtrlPts, WheelRadius, MarkerRadius, ...
  MarkerAngle0, Time0, MaxDistDelta )

% initial guess for time
DistDelta = 1/ceil(1/MaxDistDelta);
LocalTime = 0:DistDelta:1;

while true

% compute the points on the Bezier curve and the wheel that rolls over it
BezierPos  = EvalBezier( CtrlPts, LocalTime );
BezierNorm = EvalBezierNormal( CtrlPts, LocalTime, WheelRadius );
if WheelRadius > 0
  WhCtrPos = BezierPos + BezierNorm;
else
  WhCtrPos = BezierPos;
end
BezNormAngleDiff = diff( atan2(BezierNorm(2,:), BezierNorm(1,:)), 1, 2 );

% angle that the wheel spun between two given points
if WheelRadius > 0
  DiffAngle = vecnorm( diff(BezierPos,1,2), 2, 1) / WheelRadius;
else
  DiffAngle = vecnorm( diff(BezierPos,1,2), 2, 1);
end

% position and angle for the marker
MarkerAngle = cumsum([MarkerAngle0, -DiffAngle+BezNormAngleDiff]);
MarkerPos   = WhCtrPos + [cos(MarkerAngle); sin(MarkerAngle)]*MarkerRadius;

% check if the marker points are not too far from each other
DiffCurve = vecnorm( diff(MarkerPos,1,2), 2, 1);
if max(DiffCurve) < MaxDistDelta
  break
end

% if the marker points are too far, add more time points when needed
NewTimes = [];
for i = 2:length(LocalTime)
  if( DiffCurve(i-1) > MaxDistDelta )
    epsilon = (LocalTime(i)-LocalTime(i-1))/ceil(DiffCurve(i-1)/(MaxDistDelta/2));
    NewTimes = [ NewTimes, (LocalTime(i-1):epsilon:LocalTime(i)) ];
  end
end

% use the new timepoints and iterate until the result is acceptable
LocalTime = unique( [LocalTime, NewTimes], "sorted" );
end

Time = LocalTime + Time0;

% technical
MarkerAngle = mod(MarkerAngle, 2*pi);

end

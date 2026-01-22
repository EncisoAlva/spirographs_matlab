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
%           Tol  Maximum allowable distance between neighboring points [1]
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
function [Time, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  SingleBezierSegment_gral( CtrlPts, WheelRadius, MarkerRadius, ...
  MarkerAngle0, Time0, Tol, ExtraArgs )

% initial guess for time
DistDelta = 1/ceil(1/Tol);
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
switch ExtraArgs.Wheeltype
  case 'polygon'
    aang = 2*pi/ExtraArgs.nSides;
    SideAngle = mod(MarkerAngle, aang);
    %
    X0 = -sin( -SideAngle )*MarkerRadius;
    Y0 = -cos( -SideAngle )*MarkerRadius;
    X1 = -sin( -SideAngle + aang )*MarkerRadius;
    Y1 = -cos( -SideAngle + aang )*MarkerRadius;
    %
    MarkerPos = WhCtrPos + ...
      ( (Y1-Y0)*X0 - (X1-X0)*Y0 )./( (Y1-Y0)*0 - (X1-X0)*1 );
  case 'polygon_static'
    %
    MarkerAngle_sh = mod(MarkerAngle, 2*pi);
    MarkerPos = WhCtrPos*0;
    idx_tmp = 1:size(MarkerAngle,2);
    for s = 0:( ExtraArgs.nSides-1 )
      curr_pts = idx_tmp((MarkerAngle_sh >= s*2*pi/ExtraArgs.nSides)&(MarkerAngle_sh <= (s+1)*2*pi/ExtraArgs.nSides));
      %
      x0 = -sin(  s   *2*pi/ExtraArgs.nSides );
      y0 =  cos(  s   *2*pi/ExtraArgs.nSides );
      x1 = -sin( (s+1)*2*pi/ExtraArgs.nSides );
      y1 =  cos( (s+1)*2*pi/ExtraArgs.nSides );
      %
      MarkerPos(curr_pts) = WhCtrPos(curr_pts) + ...
        ( (y1-y0)*x0 - (x1-x0)*y0 )./( (y1-y0)*cos(MarkerAngle_sh(curr_pts)) - (x1-x0)*sin(MarkerAngle_sh(curr_pts)) );
    end
  otherwise
    MarkerPos   = WhCtrPos + [cos(MarkerAngle); sin(MarkerAngle)]*MarkerRadius;
end

% check if the marker points are not too far from each other
DiffCurve = vecnorm( diff(MarkerPos,1,2), 2, 1);
if max(DiffCurve) < Tol
  break
end

% if the marker points are too far, add more time points when needed
NewTimes = [];
for i = 2:length(LocalTime)
  if( DiffCurve(i-1) > Tol )
    epsilon = (LocalTime(i)-LocalTime(i-1))/ceil(DiffCurve(i-1)/(Tol/2));
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

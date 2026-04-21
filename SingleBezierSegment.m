% Create the spirograph for one Bezier cubic curve; the full curve is the
% concatenation of multiple Bezier curves.
%
% ---- INUPUT ------------------------------------------------------------
%       CtrlPts  Control points for Bezier curve, [2x4]
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                wheelcenter-marker line [1]
%         Time0  Initial timestamp [1]
%     ExtraOpts  More arguments, including optional [structs]
% -------------
%           Tol  Maximum allowable distance between neighboring points [1]
%        Method  What is inside the rolling circle
%  --->  Method = Default : one single point
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%  --->  Method = Hole
%  --->  Method = Ring
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
  SingleBezierSegment( CtrlPts, WheelRadius, MarkerAngle0, Time0, ...
  ExtraOpts )

% dealing with optional arguments
if ~isfield(ExtraOpts, 'Tol')
  Tol = 1e-3; 
else
  Tol = ExtraOpts.Tol;
end
if ~isfield(ExtraOpts, 'Method')
  ExtraOpts.Method = 'Default';
end
switch ExtraOpts.Method
  case 'Default'
    if ~isfield(ExtraOpts, 'MarkerRadius')
      MarkerRadius = 0;
    else
      MarkerRadius = ExtraOpts.MarkerRadius;
    end
  case 'Hole'
    if ~isfield(ExtraOpts, 'RollDist0')
      RollDist0 = 0; 
    else
      RollDist0 = ExtraOpts.RollDist0; 
    end
    BezBase = ExtraOpts.BezBase;
    AngBase = ExtraOpts.AngBase;
  case 'Ring2'
    if ~isfield(ExtraOpts, 'RollDist0')
      RollDist0 = 0; 
    else
      RollDist0 = ExtraOpts.RollDist0; 
    end
    CtrHoleDist   = ExtraOpts.CtrHoleDist;
    HoleRadius    = ExtraOpts.HoleRadius;
    Wheel2Radius  = ExtraOpts.Wheel2Radius;
    Marker2Radius = ExtraOpts.Marker2Radius;
end

% initial guess for time
PerUpBound = ...
  norm(CtrlPts(:,1)-CtrlPts(:,2)) + ...
  norm(CtrlPts(:,2)-CtrlPts(:,3)) + ...
  norm(CtrlPts(:,3)-CtrlPts(:,4));
DistDelta = 1/ceil(PerUpBound/Tol);
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
switch ExtraOpts.Method
  case 'Default'
    % add a point in a circle with given center and radius
    MarkerPos = WhCtrPos + [cos(MarkerAngle); sin(MarkerAngle)]*MarkerRadius;
  case 'Hole'
    % distance that the wheel has rolled so far
    RollAngle = mod( cumsum([RollDist0, DiffAngle ]), 2*pi);
    % 1. interpolate where the marker is in the hole
    HolePos = zeros(size(WhCtrPos));
    HolePos(1,:) = interp1(AngBase, BezBase(1,:), mod(RollAngle,2*pi));
    HolePos(2,:) = interp1(AngBase, BezBase(2,:), mod(RollAngle,2*pi));
    % 2. rotate the interpolated point and add around the wheel center
    MarkerPos = zeros(size(WhCtrPos));
    for j = 1:size(WhCtrPos,2)
      th = MarkerAngle(j);
      MarkerPos(:,j) = WhCtrPos(:,j) + ...
        WheelRadius * [cos(th) -sin(th); sin(th) cos(th)] * HolePos(:,j);
    end
  case 'Ring2'
    % distance that the wheel has rolled so far
    RollAngle = mod( cumsum([RollDist0, DiffAngle ]), 2*pi);
    % compute the angle rolled by the smallest gear
    CircProject = [cos(RollAngle); sin(RollAngle)] - [CtrHoleDist;0];
    LocRollAngle = atan2(CircProject(2,:),CircProject(1,:));
    % position of ring center IF the wheel was static
    LocHoleCtrPos = WhCtrPos + CtrHoleDist * [cos(MarkerAngle); sin(MarkerAngle)];
    % position of marker IF ring was static
    LocMarkerPos = [...
      (Wheel2Radius-HoleRadius)*cos(LocRollAngle) + Marker2Radius*cos(LocRollAngle*((Wheel2Radius-HoleRadius)/HoleRadius)),...
      (Wheel2Radius-HoleRadius)*sin(LocRollAngle) + Marker2Radius*sin(LocRollAngle*((Wheel2Radius-HoleRadius)/HoleRadius))...
      ];
    % rotate the interpolated point and add around the wheel center
    disp('')
    MarkerPos = zeros(size(WhCtrPos));
    for j = 1:size(WhCtrPos,2)
      th = MarkerAngle(j);
      MarkerPos(:,j) = LocHoleCtrPos(:,j) + ...
        WheelRadius * [cos(th) -sin(th); sin(th) cos(th)] * LocMarkerPos(:,j);
    end
  otherwise
    MarkerPos = WhCtrPos;
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

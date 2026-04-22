% Create the spirograph for multiple concatenated Bezier curves.
%
% ---- INUPUT ------------------------------------------------------------
%         BPath  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%   WheelRadius  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%  MarkerRadius  Distance from the center of wheel to the marker; if it is
%                larger than the wheel radius, the marker is outside [1]
%  MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                wheelcenter-marker line [1]
%      CloseTol  Max distance between first and last point [1]
%      MaxSpins  Max full rotations of wheel around whole shape [1]
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
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [Time, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette( BPath, WheelRadius, MarkerAngle0, ...
    ExtraOpts)

% reading aditional parameters
if ~isfield(ExtraOpts, 'Method')
  ExtraOpts.Method = 'Default';
end
SegmentOpts = [];
SegmentOpts.Tol = ExtraOpts.Tol;
SegmentOpts.Method = ExtraOpts.Method;
switch ExtraOpts.Method
  case 'Default'
    if ~isfield(ExtraOpts, 'MarkerRadius')
      SegmentOpts.MarkerRadius = 0;
    else
      SegmentOpts.MarkerRadius = ExtraOpts.MarkerRadius;
    end
  case 'Hole'
    SegmentOpts.BezBase = ExtraOpts.BezBase;
    SegmentOpts.AngBase = ExtraOpts.AngBase;
  case 'Ring2'
    SegmentOpts.CtrHoleDist   = ExtraOpts.CtrHoleDist;
    SegmentOpts.HoleRadius    = ExtraOpts.HoleRadius;
    SegmentOpts.Wheel2Radius  = ExtraOpts.Wheel2Radius;
    SegmentOpts.Marker2Radius = ExtraOpts.Marker2Radius;
end
%
if isfield(ExtraOpts,'MinSpins')
  MinSpins = ExtraOpts.MinSpins;
else
  MinSpins = 0;
end
if isfield(ExtraOpts,'MaxSpins')
  MaxSpins = ExtraOpts.MaxSpins;
else
  MaxSpins = 100;
end
if isfield(ExtraOpts,'Tol')
  Tol = ExtraOpts.Tol;
else
  maxX = -Inf;
  minX =  Inf;
  maxY = -Inf;
  minY =  Inf;
  for p = 1:length(BPath)
    minX = min( minX, min(BPath{p}(1,:)) );
    maxX = max( maxX, max(BPath{p}(1,:)) );
    minY = min( minY, min(BPath{p}(2,:)) );
    maxY = max( maxY, max(BPath{p}(2,:)) );
  end
  %
  Tol = norm( [ maxX-minX; maxY-minY ] )*(1e-3);
end
if isfield(ExtraOpts,'CloseTol')
  CloseTol = ExtraOpts.CloseTol;
else
  CloseTol = Tol*(1e-1);
end

% parameters
nCurves    = length(BPath);

% containers for results
Time        = [];
BezierPos   = [];
WhCtrPos    = [];
MarkerPos   = [];
MarkerAngle = [];

% initialize
FirstTangent = EvalBezierNormal(BPath{1},0,1);
CurrAngle0   = atan2(FirstTangent(2), FirstTangent(1)) + MarkerAngle0 + pi; % point to the curve
CurrTime0    = 0;

% loop
CurrSpin = 0; % index start at 0
ClosedFlag = false;
SufficientSpins = false;
CurrRollDist0 = 0;
while (CurrSpin < MaxSpins) && (~ClosedFlag)
  for j = 1:nCurves
    disp(strcat('Spin: ',num2str(CurrSpin),' , Curve: ',num2str(j)))
    CurrCtrlPts = BPath{j};
    %
    % parsing arguments
    switch ExtraOpts.Method
      case 'Default'
        % nothing
      case 'Hole'
        SegmentOpts.RollDist0 = CurrRollDist0;
      case 'Ring2'
        SegmentOpts.RollDist0 = CurrRollDist0;
    end
    % run one single Bezier curve at the time
    [locTime, locBezierPos, locWhCtrPos, locMarkerPos, locMarkerAngle] = ...
      SingleBezierSegment( CurrCtrlPts, WheelRadius, CurrAngle0, CurrTime0, ...
        SegmentOpts );
    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        locTime];
    BezierPos   = [BezierPos,   locBezierPos];
    WhCtrPos    = [WhCtrPos,    locWhCtrPos];
    MarkerPos   = [MarkerPos,   locMarkerPos];
    MarkerAngle = [MarkerAngle, locMarkerAngle];
    %
    % update initial values
    CurrTime0  = Time(end);
    CurrAngle0 = MarkerAngle(end);
    CurrRollDist0 = CurrRollDist0 + PathPerimeter(CurrCtrlPts, Tol)/WheelRadius;
    %
    % prepare for a corner
    if j<nCurves
      NextCtrlPts = BPath{j+1};
    else
      NextCtrlPts = BPath{1};
    end
    %
    % roll over the corner, if needed
    [locTime, locBezierPos, locWhCtrPos, locMarkerPos, locMarkerAngle] = ...
      RollCorner( CurrCtrlPts, NextCtrlPts, WheelRadius, MarkerPos(:,end), ...
      CurrAngle0, CurrTime0, Tol );
    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        locTime];
    BezierPos   = [BezierPos,   locBezierPos];
    WhCtrPos    = [WhCtrPos,    locWhCtrPos];
    MarkerPos   = [MarkerPos,   locMarkerPos];
    MarkerAngle = [MarkerAngle, locMarkerAngle];
    %
  %   % debug
  %   if ExtraOpts.Method == 'Ring2'
  %     LocHoleCtrPos = WhCtrPos + SegmentOpts.CtrHoleDist * [cos(MarkerAngle); sin(MarkerAngle)];
  %     figure()
  %     hold on
  %     axis equal
  %     plot(BezierPos(1,:),BezierPos(2,:))
  %     plot(WhCtrPos(1,:),WhCtrPos(2,:))
  %     plot(LocHoleCtrPos(1,:),LocHoleCtrPos(2,:))
  %     plot(MarkerPos(1,:),MarkerPos(2,:))
  %   end
  %   %
  %   % update initial values
  %   CurrTime0  = Time(end);
  %   CurrAngle0 = MarkerAngle(end);
  % end
  %
  % update number of spins
  CurrSpin = CurrSpin + 1;
  if CurrSpin >= MinSpins
    SufficientSpins = true;
  end
  %
  % check if spirograph is closed
  if CurrSpin == 1
    FirstPt = MarkerPos(:,1);
  end
  LastPt = MarkerPos(:,end);
  if norm( FirstPt - LastPt ) < CloseTol
    if SufficientSpins
      ClosedFlag = true;
    end
  end
end

end
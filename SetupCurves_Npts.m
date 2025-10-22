% Given a list of Bezier curves and parameters, creates arrays 
% necessary to create the animations.
% The goal is to test multiple parameter configurations before making the 
% animation, since testing parameters is much faster than the rest of the
% animation.
%
% ---- INUPUT ------------------------------------------------------------
%           BPath  Array with control points for each one of the Bezier
%                  curves that make the curve {?} <- [2,4]'s
%   WheelBezRatio  Perimeter of shape / Perimeter of wheel [1]
% WheelMarkerRatio Center of wheel to marker / Radius of wheel [1]
%    MarkerAngle0  Initial angle between the wheelcenter-curve line and the
%                  the wheelcenter-marker line [1]
%             Tol  Maxi allowable distance between neighboring points [1]
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
  SetupCurves_Npts( nPts,...
  BPath, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
    MaxSpins, ExtraOpts)

% this evaluation is for background decoration only
DecorativeBez = PathEval(BPath, ExtraOpts.Tol);

% containers
AllBezierPos   = cell(1,nPts);
AllLocTime     = cell(1,nPts);
AllWhCtrPos    = cell(1,nPts);
AllMarkerPos   = cell(1,nPts);
AllMarkerAngle = cell(1,nPts);

% loop to create multiple curves
for i = 1:nPts
  [LocTime, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
    GenerateGlissette( BPath, WheelRadius, MarkerRadius, MarkerAngle0Array(i), ...
    ExtraOpts.Tol, ExtraOpts.CloseTol, MaxSpins);

  % patch
  if ExtraOpts.CloseEnds
    MarkerPos(:,end+1) = MarkerPos(:,1);
    MarkerAngle(end+1) = MarkerAngle(1);
  end

  % reverse curves with opposite orientation
  if isfield(ExtraOpts,'ChangeOrient')
  if ExtraOpts.ChangeOrient(i)
    BezierPos   = flip(BezierPos, 2);
    WhCtrPos    = flip( WhCtrPos, 2 );
    MarkerPos   = flip(MarkerPos, 2);
    MarkerAngle = flip(MarkerAngle, 2);
    LocTime = max(LocTime) - flip(LocTime);
  end
  end

  % Store results in containers
  AllBezierPos{i} = BezierPos;
  AllLocTime{i} = LocTime;
  AllWhCtrPos{i} = WhCtrPos;
  AllMarkerPos{i} = MarkerPos;
  AllMarkerAngle{i} = MarkerAngle;
end

end
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
function [ DecorativeBez,...
  AllBezierPos, AllLocTime, ...
  AllWhCtrPos, AllMarkerPos, AllMarkerAngle ] = ...
  SetupCurves_4pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, ...
    CloseTol, MaxSpins)

% this evaluation is for background decoration only
DecorativeBez = AllBezierEval(CtrlPtsArray, MaxDistDelta);

% containers
AllBezierPos   = {};
AllLocTime     = {};
AllWhCtrPos    = {};
AllMarkerPos   = {};
AllMarkerAngle = {};

% outer curves
[BezierPos1A, BezierPos2A, ...
  LocTime1A,...
  WhCtrPos1A, MarkerPos1A, MarkerAngle1A,...
  LocTime2A,...
  WhCtrPos2A, MarkerPos2A, MarkerAngle2A] = ...
  SetupCurves_2pts( CtrlPtsArray, WheelRadius, MarkerRadius, MarkerAngle0, ...
    MaxDistDelta, CloseTol, MaxSpins);

[BezierPos1B, BezierPos2B, ...
  LocTime1B,...
  WhCtrPos1B, MarkerPos1B, MarkerAngle1B,...
  LocTime2B,...
  WhCtrPos2B, MarkerPos2B, MarkerAngle2B] = ...
  SetupCurves_2pts( FlipBezierAll(CtrlPtsArray), WheelRadius, MarkerRadius, -MarkerAngle0+pi, ...
    MaxDistDelta, CloseTol, MaxSpins);

BezierPos1B   = flip(BezierPos1B, 2);
WhCtrPos1B    = flip( WhCtrPos1B, 2 );
MarkerPos1B   = flip(MarkerPos1B, 2);
MarkerAngle1B = flip(MarkerAngle1B, 2);
%
BezierPos2B   = flip(BezierPos2B, 2);
WhCtrPos2B    = flip(WhCtrPos2B, 2);
MarkerPos2B   = flip(MarkerPos2B, 2);
MarkerAngle2B = flip(MarkerAngle2B, 2);

LocTime1B = max(LocTime1B) - flip(LocTime1B);
LocTime2B = max(LocTime2B) - flip(LocTime2B);

% Store results in containers
AllBezierPos{1} = BezierPos1A;
AllBezierPos{2} = BezierPos2A;
AllBezierPos{3} = BezierPos1B;
AllBezierPos{4} = BezierPos2B;
%
AllLocTime{1} = LocTime1A;
AllLocTime{2} = LocTime2A;
AllLocTime{3} = LocTime1B;
AllLocTime{4} = LocTime2B;
%
AllWhCtrPos{1} = WhCtrPos1A;
AllWhCtrPos{2} = WhCtrPos2A;
AllWhCtrPos{3} = WhCtrPos1B;
AllWhCtrPos{4} = WhCtrPos2B;
%
AllMarkerPos{1} = MarkerPos1A;
AllMarkerPos{2} = MarkerPos2A;
AllMarkerPos{3} = MarkerPos1B;
AllMarkerPos{4} = MarkerPos2B;
%
AllMarkerAngle{1} = MarkerAngle1A;
AllMarkerAngle{2} = MarkerAngle2A;
AllMarkerAngle{3} = MarkerAngle1B;
AllMarkerAngle{4} = MarkerAngle2B;

end
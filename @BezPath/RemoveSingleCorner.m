% Corner case: the wheel can't continue rolling on the current curve
% because it will collide with the next curve.
% Under such assumption, this function will find the indexes ct1, ct2 from
% the curves C1 and C2 such that the scaled normals at C1(ct1) and C2(ct2) 
% are equal --up to some tolerance. This means that, it the curve was a
% physical object, the wheel will stop there.
% Needless to say, this process depends heavily on the wheel size.
%
% ---- INUPUT ------------------------------------------------------------
%       Curve1/2  Array with control points for each one of the Bezier
%                 curves that make the curve {?} <- [2,4]'s
%            Tol  Max distance between points of the discretization [1]
%
% ---- OUTPUT ------------------------------------------------------------
% Curve1/2_short  Array with control points for each one of the Bezier
%                 curves AFTER replacing the corner portion {?} <- [2,4]'s
%           Time  Timestamps [1x?]
%       WhCtrPos  Location of the wheel center at timestamps [2x?]
%      MarkerPos  Location of marker at timepoints [2x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [SuccessFlag, CtrlPtsPrev_new, CtrlPtsPost_new, CtrlPts_roll] = ...
  RemoveSingleCorner( CtrlPtsPrev, CtrlPtsPost, WheelRadius )

% find collision time of curves, if there is one
%SegmentPrev = BezSegment( CtrlPtsPrev );
%SegmentPost = BezSegment( CtrlPtsPost );

[ct1,ct2] = ...
  BezSegment.FindCollisionTime( CtrlPtsPrev, CtrlPtsPost, WheelRadius );

% do stuff only if there is a roll
SuccessFlag = false;
CtrlPtsPrev_new = zeros(2,4);
CtrlPtsPost_new = zeros(2,4);
CtrlPts_roll    = zeros(2,4);
if ~isnan(ct1) && ~isnan(ct2)

% adjusting old format, may rewrtite later
SuccessFlag = true;

% Casteljau division at the point of intersection
CtrlPtsPrev_new(:,1) = CtrlPtsPrev(:,1);
CtrlPtsPrev_new(:,2) = (1-ct1)*CtrlPtsPrev(:,1) + ct1*CtrlPtsPrev(:,2);
CtrlPtsPrev_new(:,3) = (1-ct1)^2*CtrlPtsPrev(:,1) + 2*(1-ct1)*ct1*CtrlPtsPrev(:,2) + ct1^2*CtrlPtsPrev(:,3);
CtrlPtsPrev_new(:,4) = (1-ct1)^3*CtrlPtsPrev(:,1) + 3*(1-ct1)^2*ct1*CtrlPtsPrev(:,2) + 3*(1-ct1)*ct1^2*CtrlPtsPrev(:,3) + ct1^3*CtrlPtsPrev(:,4);

CtrlPtsPost_new(:,4) = CtrlPtsPost(:,4);
CtrlPtsPost_new(:,3) = (1-ct2)*CtrlPtsPost(:,3) + ct2*CtrlPtsPost(:,4);
CtrlPtsPost_new(:,2) = (1-ct2)^2*CtrlPtsPost(:,2) + 2*(1-ct2)*ct2*CtrlPtsPost(:,3) + ct2^2*CtrlPtsPost(:,4);
CtrlPtsPost_new(:,1) = (1-ct2)^3*CtrlPtsPost(:,1) + 3*(1-ct2)^2*ct2*CtrlPtsPost(:,2) + 3*(1-ct2)*ct2^2*CtrlPtsPost(:,3) + ct2^3*CtrlPtsPost(:,4);

% getting normal vectors to use the 4/3 tan(th/4) approx of a circle
SegmentPrev_new = BezSegment( CtrlPtsPrev_new );
SegmentPost_new = BezSegment( CtrlPtsPost_new );
N1 = SegmentPrev_new.EvalNormal( ct1, 1 );
N2 = SegmentPost_new.EvalNormal( ct2, 1 );

GapAngle = atan2(N2(2),N2(1)) - atan2(N1(2),N1(1));
if GapAngle < 0
  GapAngle = -( 2*pi - GapAngle );
end

% gap is filled with a circle arc
CtrlPts_roll(:,1) = SegmentPrev_new.EvalPosition( ct1 );
CtrlPts_roll(:,4) = SegmentPost_new.EvalPosition( ct2 );
CtrlPts_roll(:,2) = CtrlPts_roll(:,1) + [0,1;-1,0]*N1* (4/3)*tan(GapAngle/4)*WheelRadius;
CtrlPts_roll(:,3) = CtrlPts_roll(:,4) + [0,-1;1,0]*N2* (4/3)*tan(GapAngle/4)*WheelRadius;

end

end
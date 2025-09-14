% Corner case: the wheel can't continue rolling on the current curve
% because it will collide with the next curve.
% Under such assumption, this function will find the indexes t1, t2 from
% the curves C1 and C2 such that the scaled normals at C1(t1) and C2(t2) 
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
function [Curve1_short, Curve2_short, Curve_gap] = ...
  RemoveCorner( CurvePre, CurvePos, WheelRadius, Tol )

[t1,t2] = FindCollisionTime( CurvePre, CurvePos, WheelRadius, Tol );

% Casteljau division at the point of intersection
Curve1_short = zeros(2,4);
Curve2_short = zeros(2,4);

Curve1_short(:,1) = CurvePre(:,1);
Curve1_short(:,2) = (1-t1)*CurvePre(:,1) + t1*CurvePre(:,2);
Curve1_short(:,3) = (1-t1)^2*CurvePre(:,1) + 2*(1-t1)*t1*CurvePre(:,2) + t1^2*CurvePre(:,3);
Curve1_short(:,4) = (1-t1)^3*CurvePre(:,1) + 3*(1-t1)^2*t1*CurvePre(:,2) + 3*(1-t1)*t1^2*CurvePre(:,3) + t1^3*CurvePre(:,4);

Curve2_short(:,4) = CurvePos(:,4);
Curve2_short(:,3) = (1-t2)*CurvePos(:,3) + t2*CurvePos(:,4);
Curve2_short(:,2) = (1-t2)^2*CurvePos(:,2) + 2*(1-t2)*t2*CurvePos(:,3) + t2^2*CurvePos(:,4);
Curve2_short(:,1) = (1-t2)^3*CurvePos(:,1) + 3*(1-t2)^2*t2*CurvePos(:,2) + 3*(1-t2)*t2^2*CurvePos(:,3) + t2^3*CurvePos(:,4);

% getting normal vectors to use the 4/3 tan(th/4) approx of a circle
N1 = EvalBezierNormal(CurvePre,t1,1);
N2 = EvalBezierNormal(CurvePos,t2,1);
GapAngle = atan2(N2(2),N2(1)) - atan2(N1(2),N1(1));
if GapAngle < 0
  GapAngle = -( 2*pi - GapAngle );
end

% gap is filled with a circle arc
Curve_gap = zeros(2,4);

Curve_gap(:,1) = EvalBezier(CurvePre,t1);
Curve_gap(:,4) = EvalBezier(CurvePos,t2);
Curve_gap(:,2) = Curve_gap(:,1) + [0,1;-1,0]*N1* (4/3)*tan(GapAngle/4)*WheelRadius;
Curve_gap(:,3) = Curve_gap(:,4) + [0,-1;1,0]*N2* (4/3)*tan(GapAngle/4)*WheelRadius;

end
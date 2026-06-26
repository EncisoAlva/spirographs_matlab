% Corner case: the wheel can't continue rolling on the current curve
% because it will collide with the next curve.
% Under such assumption, this function will find the indexes t1, t2 from
% the curves C1 and C2 such that the scaled normals at C1(t1) and C2(t2) 
% are equal --up to some tolerance. This means that, it the curve was a
% physical object, the wheel will stop there.
% Needless to say, this process depends heavily on the wheel size.
%
% ---- INUPUT ------------------------------------------------------------
%      Curve1/2  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%           Tol  Max distance between points of the discretization [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      WhCtrPos  Location of the wheel center at timestamps [2x?]
%     MarkerPos  Location of marker at timepoints [2x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [CrossTime1, CrossTime2] = ...
  FindCollisionTime( Curve1, Curve2, WheelRadius, Tol )

% initial and finishing times
t0_1 = 0;
t0_2 = 0;
tF_1 = 1;
tF_2 = 1;

% approximate the curve: if precision up to Tol is not reached, iterate
iter = 1;
while true
  % interpolate curves adn their normals
  LocalTime1 = t0_1:((tF_1-t0_1)/(ceil((tF_1-t0_1)/Tol)+iter)):tF_1;
  LocalTime2 = t0_2:((tF_2-t0_2)/(ceil((tF_2-t0_2)/Tol)+iter)):tF_2;

  BezierPos1  = EvalBezier( Curve1, LocalTime1 );
  BezierPos2  = EvalBezier( Curve2, LocalTime2 );
  BezierNorm1 = EvalBezierNormal( Curve1, LocalTime1, WheelRadius );
  BezierNorm2 = EvalBezierNormal( Curve2, LocalTime2, WheelRadius );

  WhCtrPos1 = BezierPos1 + BezierNorm1;
  WhCtrPos2 = BezierPos2 + BezierNorm2;

  % look for the intersection per se
  curve1_dis = zeros(size(LocalTime1));
  curve1_dis_idx = zeros(size(LocalTime1));
  for i = 1:size(LocalTime1,2)
    [dist, idxx] = min( vecnorm( WhCtrPos1(:,i) - WhCtrPos2, 2, 1 ) );
    curve1_dis(i) = dist;
    curve1_dis_idx(i) = idxx;
  end

  % get the intersection time
  [dist, idx1] = min(curve1_dis);
  idx2 = curve1_dis_idx(idx1);

  % if the normals are different, make a better approximation
  if dist < Tol
    break
  end
  t0_1 = max(0,LocalTime1(idx1)-0.5/2^iter);
  t0_2 = max(0,LocalTime2(idx2)-0.5/2^iter);
  tF_1 = min(1,LocalTime1(idx1)+0.5/2^iter);
  tF_2 = min(1,LocalTime2(idx2)+0.5/2^iter);
  iter = iter+1;
end

% report results
CrossTime1 = LocalTime1(idx1);
CrossTime2 = LocalTime2(idx2);

end
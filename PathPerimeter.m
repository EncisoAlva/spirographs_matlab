% Approximate the total perimeter of a shape whose sides are cubic Bezier
% curves. The curves are approximated by a finite number of points, then
% the distance between neighboring points is computed.
%
% ---- INUPUT ------------------------------------------------------------
%         BPath  Array with control points for each one of the Bezier
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
function [Perimeter] = PathPerimeter( BPath, Tol)

Perimeter = 0;
nCurves   = length(BPath);

% loop
for j = 1:nCurves
  CurrCtrlPts = BPath{j};
  LocalTime   = 0:( 1/ceil(1/(Tol/2)) ):1;
  %iter = 0;
  while true
    BezierVals  = EvalBezier( CurrCtrlPts, LocalTime );
    DiffCurve   = vecnorm( diff(BezierVals,1,2), 2, 1);
    if max(DiffCurve) < Tol
      break
    end
    NewTimes = [];
    for i = 2:length(LocalTime)
      if( DiffCurve(i-1) > Tol )
        epsilon  = (LocalTime(i)-LocalTime(i-1))/ceil(DiffCurve(i-1)/(Tol/2));
        NewTimes = [ NewTimes, (LocalTime(i-1):epsilon:LocalTime(i)) ];
      end
    end
    LocalTime = unique( [LocalTime, NewTimes], "sorted" );
    %iter = iter +1; % additional penalization
  end
  Perimeter   = Perimeter + sum(DiffCurve);
end

end
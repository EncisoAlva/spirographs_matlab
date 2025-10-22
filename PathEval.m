% Evaluate multiple concatenated cubic Bezier curves given its control 
% points.
%
% ---- INUPUT ------------------------------------------------------------
%         BPath  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%           Tol  Max distance between points of the discretization [1]
%
% ---- OUTPUT ------------------------------------------------------------
%    BezierVals  Points in the Bezier curve [2x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [BezierVals] = PathEval( BPath, Tol)

nCurves    = length(BPath);
BezierVals = [];

% loop
for j = 1:nCurves
  CurrCtrlPts = BPath{j};
  LocalTime   = 0:( 1/ceil(1/(Tol/2)) ):1;
  %iter = 0;
  while true
    LocalBezV = EvalBezier( CurrCtrlPts, LocalTime );
    DiffCurve = vecnorm( diff(LocalBezV,1,2), 2, 1);
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
  BezierVals = [BezierVals, LocalBezV];
end

end
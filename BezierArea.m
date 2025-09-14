% Approximate the total area of a shape whose sides are cubic Bezier
% curves. The curves are approximated by a finite number of points, then
% the area of the resulting polygon is reported.
%
% ---- INUPUT ------------------------------------------------------------
%  CtrlPtsArray  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%           Tol  Max distance between points of the discretization [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Area  Area of the closed polygon used to approxiate the Bezier
%                curves [1]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [Area] = BezierArea( CtrlPtsArray, Tol)

nCurves       = length(CtrlPtsArray);
AllBezierVals = [];

% loop
for j = 1:nCurves
  CurrCtrlPts = CtrlPtsArray{j};
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
  AllBezierVals = [AllBezierVals, BezierVals];
end

Area = polyarea( AllBezierVals(1,:), AllBezierVals(2,:) );

end
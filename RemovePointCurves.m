% Change the direction of a Bezier curve.
%
% ---- INUPUT ------------------------------------------------------------
%  CtrlPtsArray  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%
% ---- OUTPUT ------------------------------------------------------------
%  CtrlPtsArrayFlipped  Array with control points for each one of the 
%                Bezier IN REVERSE ORDER {?} <- [2,4]'s
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [CtrlPtsArrayClean] = RemovePointCurves( CtrlPtsArray, Tol )

nCurves = size(CtrlPtsArray,2);
CtrlPtsArrayClean = {};

for i = 1:nCurves
  CurrCurve     = CtrlPtsArray{i};
  CurrCurveDiff = CurrCurve - mean(CurrCurve, 2);
  if max( vecnorm( CurrCurveDiff, 2, 1 ) ) > Tol
    CtrlPtsArrayClean{end+1} = CurrCurve;
  end
end

end
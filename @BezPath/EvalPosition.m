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
function [BezierVals] = EvalPosition( obj, Index, Tvals )

% just pasre the value to the appropriae segment
BezierVals = obj.Segment{Index}.EvalPosition(Tvals);

end
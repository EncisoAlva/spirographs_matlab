% Change the direction of a Bezier curve.
%
% ---- INUPUT ------------------------------------------------------------
%         BPath  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%
% ---- OUTPUT ------------------------------------------------------------
%  BPathFlipped  Array with control points for each one of the 
%                Bezier IN REVERSE ORDER {?} <- [2,4]'s
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function obj = StandardPreprocess( obj )

obj = obj.RemovePointCurves();
obj = obj.ForceCubicLines();
obj = obj.FitBox( [0;0], [2;2] );
obj = obj.PointInwards();

end
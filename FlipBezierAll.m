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
function [CtrlPtsArrayFlipped] = FlipBezierAll( CtrlPtsArray)

CtrlPtsArrayFlipped = flip( CtrlPtsArray );

for i = 1:size(CtrlPtsArray,2)
  CtrlPtsArrayFlipped{i} = flip( CtrlPtsArrayFlipped{i}, 2 );
end

end
% Traslate a path by a given vector.
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
function [BPathTraslated] = TraslatePath( BPath, xy )

BPathTraslated = BPath;
 
for i = 1:size(BPath, 2)
  BPathTraslated{i} = BPath{i} + xy;
end

end
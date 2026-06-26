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
function [BPathShifted] = ShiftPath( BPath, Shift, Halfen)

% if there is no actual shift, return the input unchanged
if (Shift == 0)&&(Halfen==false)
  BPathShifted = BPath;
  return
end

nCurves = size(BPath,2);
BPathShifted = {};

for i = 1:nCurves
  j = mod( i+Shift -1, nCurves ) + 1;
  if Halfen
    if i == 1
      [c1, c2] = HalfBezierSingle( BPath{j} );
      BPathShifted{1} = c2;
    elseif i == nCurves
      BPathShifted{end+1} = BPath{j};
      BPathShifted{end+1} = c1;
    else
      BPathShifted{end+1} = BPath{j};
    end
  else
    BPathShifted{end+1} = BPath{j};
  end
end

end
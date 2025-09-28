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
function [CtrlPtsArrayShifted] = ShiftBezierAll( CtrlPtsArray, Shift, Halfen)

% if there is no actual shift, return the input unchanged
if Shift == 0
  CtrlPtsArrayShifted = CtrlPtsArray;
  return
end

nCurves = size(CtrlPtsArray,2);
CtrlPtsArrayShifted = {};

for i = 1:nCurves
  j = mod( i+Shift -1, nCurves ) + 1;
  if Halfen
    if i == 1
      [c1, c2] = HalfBezierSingle( CtrlPtsArray{j} );
      CtrlPtsArrayShifted{1} = c2;
    elseif i == nCurves
      CtrlPtsArrayShifted{end+1} = CtrlPtsArray{j};
      CtrlPtsArrayShifted{end+1} = c1;
    else
      CtrlPtsArrayShifted{end+1} = CtrlPtsArray{j};
    end
  else
    CtrlPtsArrayShifted{end+1} = CtrlPtsArray{j};
  end
end

end
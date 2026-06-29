% Change the direction of a Bezier curve.
%
% ---- INUPUT ------------------------------------------------------------
%  BPath  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%
% ---- OUTPUT ------------------------------------------------------------
%  BPathFlipped  Array with control points for each one of the 
%                Bezier IN REVERSE ORDER {?} <- [2,4]'s
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function ForceCubicLine( obj )

% case: P0 = P1
CurrCurve = obj.CtrlPts;
if norm(CurrCurve(:,1)-CurrCurve(:,2)) < obj.Tol
  if abs( norm(CurrCurve(:,2)-CurrCurve(:,3)) + norm(CurrCurve(:,3)-CurrCurve(:,4)) - norm(CurrCurve(:,2)-CurrCurve(:,4)) ) < obj.Tol
    obj.CtrlPts = LineToBezier( CurrCurve(:,1), CurrCurve(:,4) );
  end
end

% case: P2 = P3
CurrCurve = obj.CtrlPts;
if norm(CurrCurve(:,3)-CurrCurve(:,4)) < obj.Tol
  if abs( norm(CurrCurve(:,1)-CurrCurve(:,2)) + norm(CurrCurve(:,2)-CurrCurve(:,3)) - norm(CurrCurve(:,1)-CurrCurve(:,3)) ) < obj.Tol
    obj.CtrlPts = LineToBezier( CurrCurve(:,1), CurrCurve(:,4) );
  end
end

end
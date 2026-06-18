% Evaluate a cubic Bezier curve given its control points and a set of
% values between 0 and 1, inclusive.
%
% ---- INUPUT ------------------------------------------------------------
%       CtrlPts  Control points for Bezier curve, [2x4]
%         Tvals  Values between 0 an 1, inclusive [1x?]
%
% ---- OUTPUT ------------------------------------------------------------
%    BezierVals  Points in the Bezier curve [2x?]
%
function [BezierVals] = EvalPosition( obj, TVals)
  BezierVals = ( (1-TVals).^3) .* obj.CtrlPts(:,1) + ...
    3 * (( (1-TVals).^2).*(TVals)) .* obj.CtrlPts(:,2) + ...
    3 * (( (1-TVals)).*(TVals.^2)) .* obj.CtrlPts(:,3) + ...
    ( (TVals).^3) .* obj.CtrlPts(:,4) ;
end
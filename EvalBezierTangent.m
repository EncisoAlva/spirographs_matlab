% Compute the tangent vector to a cubic Bezier curve at a set of values
% between 0 and 1.
% In order to avoid redundant computations, the points on the Bezier curve
% are taken as arguments (computed on a previous step).
%
% ---- INUPUT ------------------------------------------------------------
%       CtrlPts  Control points for Bezier curve, [2x4]
%         Tvals  Values between 0 an 1, inclusive [1x?]
%    BezierVals  Points in the Bezier curve [2x?]
%     CirRadius  Location of the wheel center at values Tvals [2x?]
%
% ---- OUTPUT ------------------------------------------------------------
%   TangentVals  Radius of spirograph wheel, a negative radius indicates
%                that the wheel rolls inside the curve [1]
%
function [BezierNormal] = EvalBezierTangent( CtrlPts, TVals, BezierVals, CirRadius)
  % the tangent vector is the derivative of the curve
  BezierPrime = 3*( (1-TVals).^2) .* (CtrlPts(:,2)-CtrlPts(:,1))+ ...
    6 * ((1-TVals).*(TVals)) .* (CtrlPts(:,3)-CtrlPts(:,2)) + ...
    3 * ( (TVals).^2) .* (CtrlPts(:,4)-CtrlPts(:,3)) ;
  % the normal tangent vector is obtained by rotating the tangent vector
  BezierNormal = CirRadius * [ -BezierPrime(2,:); BezierPrime(1,:)] ./ vecnorm(BezierPrime,2,1);
end
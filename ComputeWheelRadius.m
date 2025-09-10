% For simplicity, the ratio between the Bezier curve perimeter and the
% wheel radius must be a rational number. The properties of the resulting
% curve depends drectly from such ratio.
% Also for simplicity, the Bezier curve perimeter can be considered as
% given (especially for complicated figures) and the ratio will be adjusted
% for artistic purposes.
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
function [BezierNormal] = ComputeWheelRadius( CtrlPts, TVals, CirRadius)
  % the tangent vector is the derivative of the curve
  BezierPrime = 3*( (1-TVals).^2) .* (CtrlPts(:,2)-CtrlPts(:,1))+ ...
    6 * ((1-TVals).*(TVals)) .* (CtrlPts(:,3)-CtrlPts(:,2)) + ...
    3 * ( (TVals).^2) .* (CtrlPts(:,4)-CtrlPts(:,3)) ;
  % the normal tangent vector is obtained by rotating the tangent vector
  BezierNormal = CirRadius * [ -BezierPrime(2,:); BezierPrime(1,:)] ./ vecnorm(BezierPrime,2,1);
end
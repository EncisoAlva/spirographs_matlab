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
function [BezierNormal] = EvalBezierNormal( CtrlPts, TVals, CirRadius)
  % the tangent vector is the derivative of the curve
  BezierPrime = 3*( (1-TVals).^2) .* (CtrlPts(:,2)-CtrlPts(:,1))+ ...
    6 * ((1-TVals).*(TVals)) .* (CtrlPts(:,3)-CtrlPts(:,2)) + ...
    3 * ( (TVals).^2) .* (CtrlPts(:,4)-CtrlPts(:,3)) ;
  % the normal tangent vector is obtained by rotating the tangent vector
  BezierNormal = CirRadius * [0,-1; 1,0]*BezierPrime ./ vecnorm(BezierPrime,2,1);
  
  % interpolation for vanishing derivatives
  if any(vecnorm(BezierPrime,2,1) < 0.001)
    for i = 1:size(TVals,2)
      if norm(BezierPrime(:,i)) < 0.001
        TVals2 = ((-0.001):0.0001:0.001) + TVals(i);
        TVals2 = TVals2( (TVals2~=TVals(i)) & (TVals2>0) & (TVals2<1) );
        BezierPrime2 = 3*( (1-TVals2).^2) .* (CtrlPts(:,2)-CtrlPts(:,1))+ ...
          6 * ((1-TVals2).*(TVals2)) .* (CtrlPts(:,3)-CtrlPts(:,2)) + ...
          3 * ( (TVals2).^2) .* (CtrlPts(:,4)-CtrlPts(:,3)) ;
        BezierNormal2 = CirRadius * [0,-1; 1,0]*BezierPrime2 ./ vecnorm(BezierPrime2,2,1);
        weights = exp( -( ( TVals2-TVals(i) ).^2)/(2*(0.0005)^2) );
        BezierNormal_replacement = BezierNormal2 * weights' / sum(weights);
        BezierNormal(:,i) = BezierNormal_replacement / norm(BezierNormal_replacement);
      end
    end
  end
end
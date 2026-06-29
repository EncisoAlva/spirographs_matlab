% Take all Bezier curves and divide each one in half. The purpose is to
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
function [CtrlPts1, CtrlPts2] = HalfBezierSingle( CtrlPts)

% original
P = CtrlPts;
% first and second
Q = zeros(2,4);
R = zeros(2,4);

%  Casteljau algorithm for r = 1/2
Q(:,1) = P(:,1);
R(:,4) = P(:,4);

Q(:,2) = ( P(:,1) + P(:,2) )/2;
delta  = ( P(:,2) + P(:,3) )/2;
R(:,3) = ( P(:,3) + P(:,4) )/2;

Q(:,3) = ( Q(:,2) + delta  )/2;
R(:,2) = ( R(:,3) + delta  )/2;

Q(:,4) = ( Q(:,3) + R(:,2) )/2;
R(:,1) = Q(:,4);

CtrlPts1 = Q;
CtrlPts2 = R;

end
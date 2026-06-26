% Given 2 points, construct the line between them as a cubic Bezier curve.
%
function [CtrlPts] = LineToBezier( Pt0, PtF )

CtrlPts = zeros(2,4);

CtrlPts(:,1) = Pt0;
CtrlPts(:,2) = (2/3) * Pt0 + (1/3) * PtF;
CtrlPts(:,3) = (1/3) * Pt0 + (2/3) * PtF;
CtrlPts(:,4) = PtF;

end
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
function PointInwards( obj )

% construct convex hull
[BezPts, Curve, Tval] = obj.EvalAllPositionsExtra( 0.01 );
ConvHullIdx = convhull( BezPts' )';
ConvHullPt  = BezPts( :, ConvHullIdx );
ConvHulSegm = Curve( ConvHullIdx );
ConvHulTval = Tval( ConvHullIdx );

% check if the normal vector goes outside or inside the convex hull
CurrSegment = obj.Segment{ ConvHulSegm(1) };
NormVec = CurrSegment.EvalNormal( ConvHulTval(1), 0.1 ) + ConvHullPt(:,1);

% flip if normal vector is outside the convex hull
if ~inpolygon( NormVec(1),NormVec(2), ConvHullPt(1),ConvHullPt(2) )
  obj.Flip();
end

end
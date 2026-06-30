% Approximate the total perimeter of a shape whose sides are cubic Bezier
% curves. The curves are approximated by a finite number of points, then
% the distance between neighboring points is computed.
%
% ---- INUPUT ------------------------------------------------------------
%         BPath  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%           Tol  Max distance between points of the discretization [1]
%
% ---- OUTPUT ------------------------------------------------------------
%          Time  Timestamps [1x?]
%      WhCtrPos  Location of the wheel center at timestamps [2x?]
%     MarkerPos  Location of marker at timepoints [2x?]
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function [Perimeter] = GetPerimeter( obj )

if ~isprop( obj, 'Perimeter' ) || isempty(obj.Perimeter)

obj.Perimeter = 0;

for j = 1:obj.nSegments
  obj.Perimeter = obj.Perimeter + obj.Segment{j}.GetSegmentPerimeter();
end

end

Perimeter = obj.Perimeter;

end
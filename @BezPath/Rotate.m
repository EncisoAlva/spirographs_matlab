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
function Rotate( obj, Angle, varargin )

% by default, the center is the origin
if ~isempty(varargin)
  Center = varargin{1};
else
  Center = [0,0]';
end

for i = 1:obj.nSegments
  obj.Segment{i}.Rotate( Angle, Center );
end

end
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
function obj = Rotate(  obj, Angle, varargin )

% by default, the center of scaling is the origin
if size(varargin)>0
  Center = varargin{1};
else
  Center = [0,0]';
end

ROT = [cos(Angle), -sin(Angle); sin(Angle), cos(Angle)];
obj.CtrlPts = ROT * ( obj.CtrlPts - Center ) + Center;

end
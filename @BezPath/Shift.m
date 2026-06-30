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
function Shift( obj, Shift, varargin)

% by default, the center is the origin
if ~isempty(varargin)
  Halfen = varargin{1};
else
  Halfen = false;
end

% if there is no actual shift, return the input unchanged
if ( Shift>0 )||( Halfen==true )

if Halfen
  Segment_new = cell(1, obj.nSegments +1);
else
  Segment_new = cell(1, obj.nSegments);
end

for i = 1:obj.nSegments
  j = mod( i+Shift -1, nSegments ) + 1;
  if Halfen
    if i == 1
      [c1, c2] = obj.Segment{j}.HalfSegment();
      Segment_new{i} = BezSegment( c2 );
    else
      Segment_new{i} = obj.Segment{j};
      if i == obj.nSegments
        Segment_new{i+1} = BezSegment( c1 );
      end
    end
  else
    Segment_new{i} = obj.Segment{j};
  end
end

if Halfen
  obj.nSegments = obj.nSegments + 1;
end

end

end
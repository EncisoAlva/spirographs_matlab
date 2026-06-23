function SetColor( obj, ColorVector, varargin )
  
  if ~iscell(ColorVector)
    % single color
    obj.ColorVector = { ColorVector };
    obj.nColors = 1;
  else
    obj.ColorVector = ColorVector;
    obj.ColorRefCurve = varargin{1};
    if size(varargin)>1
      obj.ColorCycles = varargin{2};
    else
      obj.ColorCycles = 1;
    end
    obj.Multicolor = true;
    obj.nColors = size(ColorVector,2);
  end

end
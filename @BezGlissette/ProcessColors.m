% declaring common colors
function ProcessColors( obj )

if obj.Multicolor
  switch obj.ColorRefCurve
    case 'CumDist'
      ColorVal = cumsum([0, vecnorm( diff( obj.MarkerPos, 1,2), 2,1 )]);
    case 'Bezier'
      ColorVal = cumsum([0, vecnorm( diff( obj.BezierPos, 1,2), 2,1 )]);
  end
  ColorVal = ColorVal/ColorVal(end);
  ColorNum = 0.5 - 0.5*cos( obj.ColorCycles* ColorVal * 2*pi );
  %
  colorTable = zeros(obj.nColors, 3+1);
  for q = 1:obj.nColors
    colorTable(q,2:end) = obj.ColorVector;
  end
  colorTable(:,1) = linspace(0,1, obj.nColors);
  %
  obj.ColorFunc = zeros(3, size(ColorNum,2));
  obj.ColorFunc(1,:) = interp1(colorTable(:,1),colorTable(:,2), ColorNum);
  obj.ColorFunc(2,:) = interp1(colorTable(:,1),colorTable(:,3), ColorNum);
  obj.ColorFunc(3,:) = interp1(colorTable(:,1),colorTable(:,4), ColorNum);
else
  %
end

end
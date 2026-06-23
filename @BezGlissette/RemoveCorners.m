function obj = RemoveCorners( obj )
  
  % rolling both inside and outside
  if obj.RemoveCorners_Both
    WheelRadius_old = Inf;
    WheelRadius_new = (obj.OG_BPath.GetPerimeter()/(2*pi))/obj.Wheel1BezRatio;
    iter = 1;
    while (abs( WheelRadius_new - WheelRadius_old ) > obj.WheelRadiusTol) && (iter<obj.MaxIter)
      BPath_tmp = obj.OG_BPath.Flip().RemoveAllCorners( WheelRadius_new ).Flip();
      BPath_tmp = BPath_tmp.RemoveAllCorners( WheelRadius_new );
      %
      WheelRadius_old = WheelRadius_new;
      WheelRadius_new = (BPath_tmp.GetPerimeter()/(2*pi))/obj.Wheel1BezRatio;
      %
      iter = iter + 1;
    end
    obj.BPath = BPath_tmp;
  else
    % remove corners where the wheel IS NOT rolling
    if obj.RemoveCorners_NonRolling
      obj.BPath = obj.OG_BPath.Flip().RemoveAllCorners( obj.OG_BPath.CornerRoundingRadius ).Flip();
    end
    
    % remove corners where the wheel IS rolling
    if obj.RemoveCorners_Rolling
      WheelRadius_old = Inf;
      WheelRadius_new = (obj.OG_BPath.GetPerimeter()/(2*pi))/obj.Wheel1BezRatio;
      iter = 1;
      while (abs( WheelRadius_new - WheelRadius_old ) > obj.WheelRadiusTol) && (iter<obj.MaxIter)
        BPath_tmp = obj.OG_BPath.RemoveAllCorners( WheelRadius_new );
        %
        WheelRadius_old = WheelRadius_new;
        WheelRadius_new = (BPath_tmp.GetPerimeter()/(2*pi))/obj.Wheel1BezRatio;
        iter = iter + 1;
      end
      obj.BPath = BPath_tmp;
    else
      obj.BPath = obj.OG_BPath;
    end
  end

  % the 'halfen' occurs after removing corners
  obj.BPath.Shift( obj.Shiften, obj.Halfen );

  % compute wheel radius after removing corners
  obj.Wheel1Radius = (obj.BPath.GetPerimeter()/(2*pi))/obj.Wheel1BezRatio;

  % only compute marker radius if it required, it may or may not be used
  % when using a wheel with a hole, there is no marker radius
  if isprop( obj, 'WheelMarkerRatio' )
    obj.MarkerRadius = obj.Wheel1Radius*obj.WheelMarkerRatio;
  end

end
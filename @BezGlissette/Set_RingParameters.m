function obj = Set_RingParameters( obj, container )

  % reduce fraction, if needed
  G = gcd( container.RollingCurve, container.Wheel1_Outer );
  G = gcd( G, container.Wheel1_Inner );
  G = gcd( G, container.Wheel2       );
  container.RollingCurve = container.RollingCurve/G;
  container.Wheel1_Outer = container.Wheel1_Outer/G;
  container.Wheel1_Inner = container.Wheel1_Inner/G;
  container.Wheel2       = container.Wheel2      /G;

  % the ratio between the rolling curve and the wheel is quite particular
  % due to how it is handled in the rest of the code
  if obj.AutoUpdate
    obj.Set_Wheel1BezRatio( container.RollingCurve, container.Wheel1_Outer );
  else
    obj.Wheel1BezRatio_num = container.RollingCurve;
    obj.Wheel1BezRatio_den = container.Wheel1_Outer;
    obj.Wheel1BezRatio = obj.Wheel1BezRatio_num / obj.Wheel1BezRatio_den;
  end

  % catch error with zero
  if obj.Wheel1BezRatio_num * obj.Wheel1BezRatio_den == 0
    warning('Ratio of 0 or 1/0 are not allowed.')
    obj.Wheel1BezRatio_num = nan;
    obj.Wheel1BezRatio_den = nan;
    return
  end
  
  % the other parameters can be recovered from the radius of the rolling
  % wheel
  obj.Wheel1Radius_Outer = obj.Wheel1Radius;
  obj.Wheel1Radius_Inner = obj.Wheel1Radius * container.Wheel1_Inner / container.Wheel1_Outer;
  obj.Wheel2Radius = obj.Wheel1Radius * container.Wheel2 / container.Wheel1_Outer;
  obj.MarkerRadius = obj.Wheel2Radius * container.Wheel2MarkerRatio;

  if isfield(container, 'CtrHoleDist')
    % if it is given as parameter, it is given as a percent
    obj.CtrHoleDist = container.CtrHoleDist * obj.Wheel1Radius;
  else
    obj.CtrHoleDist = obj.Wheel1Radius_Outer - obj.Wheel1Radius_Inner;
  end

  % TODO: find an actual estimate for this quantity
  obj.MinSpins = 1;
  
end
function obj = Set_Wheel1BezRatio( obj, varargin )
  % read the ratio as a fraction
  if size(varargin,2) == 1
    % if only one number is given, assume it is an integer
    obj.Wheel1BezRatio_num = round( varargin{1} );
    obj.Wheel1BezRatio_den = 1;
  else
    obj.Wheel1BezRatio_num = round( varargin{1} );
    obj.Wheel1BezRatio_den = round( varargin{2} );
  end

  % error, something with zero
  if obj.Wheel1BezRatio_num * obj.Wheel1BezRatio_den == 0
    warning('Ratio of 0 or 1/0 are not allowed.')
    obj.Wheel1BezRatio_num = nan;
    obj.Wheel1BezRatio_den = nan;
    return
  end

  % reduce fraction, if needed
  G = gcd( obj.Wheel1BezRatio_num, obj.Wheel1BezRatio_den );
  obj.Wheel1BezRatio_num = obj.Wheel1BezRatio_num/G;
  obj.Wheel1BezRatio_den = obj.Wheel1BezRatio_den/G;

  % others
  obj.MinSpins = obj.Wheel1BezRatio_den;

  % for ease of use, keep both the number and the numerator and denominator
  obj.Wheel1BezRatio = obj.Wheel1BezRatio_num / obj.Wheel1BezRatio_den;
end
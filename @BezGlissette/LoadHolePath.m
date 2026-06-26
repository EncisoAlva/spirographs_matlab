function obj = LoadHolePath( obj, varargin )
  obj.OG_HPath = BezPath( varargin );
  obj.HPath = obj.OG_HPath;
end
function obj = LoadHolePath( obj, varargin )
  obj.OG_HPath = BezPath( varargin );
  obj.OG_HPath.StandardPreprocess();
  obj.HPath = obj.OG_HPath;

  % for visualization purposes
  obj.DecorativeHole = obj.HPath.EvalAllPositions();
end
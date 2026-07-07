function obj = LoadHolePath( obj, varargin )
  obj.OG_HPath = BezPath( varargin{1}, varargin{2:end} );
  obj.OG_HPath.StandardPreprocess();
  obj.OG_HPath.Flip();
  obj.HPath = obj.OG_HPath;

  % for visualization purposes
  obj.DecorativeHole = obj.HPath.EvalAllPositions();
end
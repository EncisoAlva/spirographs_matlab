function obj = LoadRollingPath( obj, varargin )
  obj.OG_BPath = BezPath( varargin );
  %obj.OG_BPath = obj.OG_BPath.Shift( obj.Shiften, obj.Halfen );
  obj.BPath = obj.OG_BPath;
end
function LoadRollingPath( obj, varargin )
  obj.OG_BPath = BezPath( varargin{1}, varargin{2:end} );
  obj.OG_BPath.StandardPreprocess();
  obj.BPath = obj.OG_BPath;
  %
  % relative tolerance
  tmpBez = obj.BPath.EvalAllPositions();
  minX = min(tmpBez(1,:));
  maxX = max(tmpBez(1,:));
  minY = min(tmpBez(2,:));
  maxY = max(tmpBez(2,:));
  obj.Tol = max( norm( [ maxX-minX; maxY-minY ] )*(1e-4), 0.0001);
  obj.CloseTol = max( obj.Tol*(1e-1), 0.0001);

  % for visualization purposes
  obj.DecorativeBez = obj.BPath.EvalAllPositions();
end
classdef BezPath
  properties
    Segment  % array, it is singular for ease of reading
    nSegments
    %
    Tol
    SkipNegPi
    %
    Area
    Perimeter
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezPath( CtrlPtsArray )
      if size( CtrlPtsArray, 2 ) > 0
        obj.nSegments = size( CtrlPtsArray, 2 );
        obj.Segment = cell( 1, obj.nSegments );
        for i = 1:obj.nSegments
          obj.Segment{i} = BezSegment( CtrlPtsArray{i} );
        end
      else
        %
      end
      obj.Tol = 0.001;
      obj.SkipNegPi = true;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = Flip( obj )
    obj = RemovePointCurves( obj )
    obj = PointInwards( obj )
    obj = FitBox( obj, Center, MaxRange )
    obj = Rescale( obj, Center, ScaleFactor )
    obj = Rotate( obj, Center, Angle )
    obj = Translate(  obj, Translation )
    %
    PlotPath( obj )
    PlotHole( obj )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Area] = GetArea( obj )
    [Perimeter] = GetPerimeter( obj )
    [BezierVals] = EvalPosition( obj )
    %
    [SuccessFlag, CtrlPtsPrev_new, CtrlPtsPost_new, CtrlPts_roll] = ...
      RemoveSingleCorner( CtrlPtsPrev, CtrlPtsPost, WheelRadius )
    obj_return = RemoveAllCorners( obj, WheelRadius )
    obj_return = RemoveAllInnerCorners( obj, WheelRadius )
    obj_return = RemoveAllOuterCorners( obj, WheelRadius )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end




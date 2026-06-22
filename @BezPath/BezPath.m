classdef BezPath
  properties
    Segment  % array, it is singular for ease of reading
    nSegments
    %
    Tol
    SkipNegPi
    CornerRoundingRadius
    %
    Area
    Perimeter
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezPath( BezPathSource, varargin )
      switch BezPathSource
        case 'CtrlPts'
          CtrlPtsArray = varargin{1};
        case 'UniqueCurve'
          CurveName = varargin{1};
          BPath_pack = struct2cell(load('ExampleCurves.mat',CurveName));
          CtrlPtsArray = BPath_pack{1};
          clear BPath_pack
        case 'IndexedCurve'
          CurveName  = varargin{1};
          CurveIndex = varargin{2};
          BPath_pack1 = struct2cell(load('ExampleCollections.mat',CurveName));
          BPath_pack2 = BPath_pack1{CurveIndex};
          CtrlPtsArray = BPath_pack2{4};
          clear BPath_pack1 BPath_pack2
      end
      %
      obj.nSegments = size( CtrlPtsArray, 2 );
      obj.Segment = cell( 1, obj.nSegments );
      for i = 1:obj.nSegments
        obj.Segment{i} = BezSegment( CtrlPtsArray{i} );
      end
      %
      obj.Tol = 0.001;
      obj.SkipNegPi = true;
      obj.Area = obj.GetArea();
      %
      obj.CornerRoundingRadius = sqrt(0.001*obj.Area/pi);
      obj.StandardPreprocess()
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = StandardPreprocess( obj )
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
    [BezierVals] = EvalPosition( obj, Tol2 )
    [BezierVals, Curve, Tval] = EvalPositionExtra( obj, Tol2 )
    %
    [SuccessFlag, CtrlPtsPrev_new, CtrlPtsPost_new, CtrlPts_roll] = ...
      RemoveSingleCorner( CtrlPtsPrev, CtrlPtsPost, WheelRadius )
    obj_return = RemoveAllCorners( obj, WheelRadius )
    obj_return = RemoveAllInnerCorners( obj, WheelRadius )
    obj_return = RemoveAllOuterCorners( obj, WheelRadius )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
    [BPathArray] = LoadSVG( filename )
    %
    function CheckExamples()
      disp('Curves included as unique curves:')
      who -file ExampleCurves.mat
      disp('Curves included as indexed curves:')
      who -file ExampleCollections.mat
    end
  end
end




classdef BezSegment
  properties
    CtrlPts
    %
    Tol
    MaxIter
    %
    Perimeter
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezSegment( CtrlPts_in )
      % 4 points  ->  curve w/ control points given
      % 2 points  ->  line, to be converted
      if size( CtrlPts_in, 2 ) == 4
        obj.CtrlPts = CtrlPts_in;
      elseif size( CtrlPts_in, 2 ) == 2
        obj.CtrlPts = zeros(2,4);
        obj.CtrlPts(:,1) = CtrlPts_in(:,1);
        obj.CtrlPts(:,2) = CtrlPts_in * [ 2/3; 1/3 ];
        obj.CtrlPts(:,3) = CtrlPts_in * [ 1/3; 2/3 ];
        obj.CtrlPts(:,4) = CtrlPts_in(:,2);
      end
      obj.Tol = 0.0001;
      obj.MaxIter = 10;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = ForceCubicLine( obj )
    obj = Flip( obj )
    obj = Rescale( obj, Center, ScaleFactor )
    obj = Rotate(  obj, Center, Angle )
    obj = Translate(  obj, Translation )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [BezierVals]    = EvalPosition( obj, TVals)
    [BezierTangent] = EvalTangent(  obj, TVals, CirRadius)
    [BezierNormal]  = EvalNormal(   obj, TVals, CirRadius)
    [CrossTime1, CrossTime2] = FindCollisionTime( obj, obj2, WheelRadius )
    [Perimeter] = GetSegmentPerimeter( obj )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end




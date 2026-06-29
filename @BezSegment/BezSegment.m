classdef BezSegment < handle
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
    ForceCubicLine( obj )
    Flip( obj )
    Rescale( obj, ScaleFactor, varargin )
    Rotate(  obj, Angle, varargin )
    Translate(  obj, Translation )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [BezierVals]    = EvalPosition( obj, TVals)
    [BezierTangent] = EvalTangent(  obj, TVals, varargin)
    [BezierNormal]  = EvalNormal(   obj, TVals, varargin)
    [Perimeter] = GetSegmentPerimeter( obj )
    [CtrlPts1, CtrlPts2] = HalfSegment( obj )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
    [CrossTime1, CrossTime2] = FindCollisionTime( CtrlPts1, CtrlPts2, WheelRadius )
  end
end




classdef BezGlissette < handle
  properties
    BPath
    OG_BPath
    Method
    %
    Wheel1BezRatio_num
    Wheel1BezRatio_den
    Wheel1BezRatio
    %
    Marker1Angle0
    Shiften
    Halfen
    RemoveCorners_Rolling
    RemoveCorners_NonRolling
    RemoveCorners_Both
    %
    ColorVector
    ColorCycles
    ColorRefCurve
    nColors
    ColorFunc
    Multicolor
    %
    Tol
    CloseTol
    MaxSpins
    WheelRadiusTol
    CloseEnds
    MaxIter
    %
    MinSpins
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezGlissette( Method )
      obj.Method = Method;
      obj.Tol      = 0.005;
      obj.CloseTol = 0.005;
      obj.WheelRadiusTol = 0.0005;
      obj.MaxSpins = 100;
      obj.MinSpins = 1;
      obj.MaxIter = 20;
      obj.Marker1Angle0 = 0;
      obj.Shiften = 0;
      obj.Halfen = false;
      obj.CloseEnds = true;
      obj.RemoveCorners_Rolling = false;
      obj.RemoveCorners_NonRolling = false;
      obj.RemoveCorners_Both = false;
      obj.Multicolor = false;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadBezPath(obj, varargin)
    Set_Wheel1BezRatio(obj, varargin)
    RemoveCorners( obj )
    SetColor( obj, ColorVector, varargin )
    ProcessColors( obj )
    PlotGlisette( obj )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
  end
end




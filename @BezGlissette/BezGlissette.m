classdef BezGlissette < handle
  properties
    Method
    BPath
    OG_BPath
    HPath
    OG_HPath
    %
    Wheel1BezRatio_num
    Wheel1BezRatio_den
    Wheel1BezRatio
    Wheel1MarkerRatio
    %
    Wheel1Radius
    MarkerRadius
    %
    MarkerAngle0
    Shiften
    Halfen
    RemoveCorners_Rolling
    RemoveCorners_NonRolling
    RemoveCorners_Both
    AutoUpdate
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
    ChangeOrient
    %
    MinSpins
    %
    DecorativeBez
    DecorativeHole
    %
    BezierPos
    WhCtrPos
    MarkerPos
    MarkerAngle
    LocTime
    %
    BezBase
    AngBase
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
      obj.MarkerAngle0 = 0;
      obj.Shiften = 0;
      obj.Halfen = false;
      obj.CloseEnds = true;
      obj.RemoveCorners_Rolling = false;
      obj.RemoveCorners_NonRolling = false;
      obj.RemoveCorners_Both = false;
      obj.Multicolor = false;
      obj.ChangeOrient = false;
      obj.AutoUpdate = false;
      %
      obj.Tol = 0.0001;
      obj.CloseTol = 0.0001;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadRollingPath(obj, varargin)
    LoadHolePath(obj, varargin)
    Set_Wheel1BezRatio(obj, varargin)
    RemoveCorners( obj )
    SetColor( obj, ColorVector, varargin )
    SetupHole( obj, varargin )
    DEV_SetupHole_concave( obj, varargin )
    ProcessColors( obj )
    PlotGlisette( obj )
    PlotStarFilling( obj )
    ProcessGlissette( obj )
    MakeVideo( obj, VidName, VideoOpts, varargin )
    %
    PlotPath( obj )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
  end
end




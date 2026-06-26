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
    %
    Wheel1Radius
    MarkerRadius
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
    ChangeOrient
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
      obj.ChangeOrient = false;
      %
      % relative tolerance
      minX = min(obj.BPath.CtrlPts(1,:));
      maxX = max(obj.BPath.CtrlPts(1,:));
      minY = min(obj.BPath.CtrlPts(2,:));
      maxY = max(obj.BPath.CtrlPts(2,:));
      obj.Tol = max( norm( [ maxX-minX; maxY-minY ] )*(1e-4), 0.005);
      obj.CloseTol = max( obj.Tol*(1e-1), 0.001);
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadRolingPath(obj, varargin)
    LoadHolePath(obj, varargin)
    Set_Wheel1BezRatio(obj, varargin)
    RemoveCorners( obj )
    SetColor( obj, ColorVector, varargin )
    SetupHole( obj, CheckFit )
    DEV_SetupHole_concave( obj, CheckFit )
    ProcessColors( obj )
    PlotGlisette( obj )
    ProcessGlissette( obj )
    MakeVideo( obj )
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
  end
end




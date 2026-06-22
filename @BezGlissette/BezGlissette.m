classdef BezGlissette < handle
  properties
    BPath
    %
    Marker1Angle0
    Shiften
    Halfen
    %
    Tol
    CloseTol
    MaxSpins
    WheelRadiusTol
    %
    MinSpins
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezGlissette()
      obj.Tol      = 0.005;
      obj.CloseTol = 0.005;
      obj.WheelRadiusTol = 0.0005;
      obj.MaxSpins = 100;
      obj.MinSpins = 1;
      obj.Marker1Angle0 = 0;
      obj.Shiften = 0;
      obj.Halfen = false;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadBezPath(obj, varargin)
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
  end
end




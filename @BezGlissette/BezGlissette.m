classdef BezGlissette < handle
  properties
    BPath
    %
    Wheel1BezRatio_num
    Wheel1BezRatio_den
    %
    Marker1Angle0
    Shiften
    Halfen
    RemoveCorners_Rolling
    RemoveCorners_NonRolling
    %
    Tol
    CloseTol
    MaxSpins
    WheelRadiusTol
    CloseEnds
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
      obj.CloseEnds = true;
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LoadBezPath(obj, varargin)
    Set_Wheel1BezRatio(obj, varargin)
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  methods (Static)
  end
end




classdef BezPath
  properties
    Segment  % array, it is singular for ease of reading
    nSegments
    %
    Tol
    %
    Area
  end

  methods
    %%%  CONSTRUCTORS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = BezPath()
      %
    end
    %%%  METHODS ; OUTPUT = NO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  METHODS ; OUTPUT = YES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Area] = GetArea( obj )
    [Perimeter] = GetPerimeter( obj )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end




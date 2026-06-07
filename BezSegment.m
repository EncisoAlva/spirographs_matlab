classdef BezSegment
  properties
    CtrlPts
    Tol
  end

  methods
    function obj = BezSegment(CtrlPts_in)
      obj.CtrlPts = CtrlPts_in;
      obj.Tol = 0.0001;
    end
    %%
    % Evaluate a cubic Bezier curve given its control points and a set of
    % values between 0 and 1, inclusive.
    %
    % ---- INUPUT ------------------------------------------------------------
    %       CtrlPts  Control points for Bezier curve, [2x4]
    %         Tvals  Values between 0 an 1, inclusive [1x?]
    %
    % ---- OUTPUT ------------------------------------------------------------
    %    BezierVals  Points in the Bezier curve [2x?]
    %
    function [BezierVals] = EvalPosition( obj, TVals)
      BezierVals = ( (1-TVals).^3) .* obj.CtrlPts(:,1) + ...
      3 * (( (1-TVals).^2).*(TVals)) .* obj.CtrlPts(:,2) + ...
      3 * (( (1-TVals)).*(TVals.^2)) .* obj.CtrlPts(:,3) + ...
      ( (TVals).^3) .* obj.CtrlPts(:,4) ;
    end
    % Compute the tangent vector to a cubic Bezier curve at a set of values
    % between 0 and 1.
    % In order to avoid redundant computations, the points on the Bezier curve
    % are taken as arguments (computed on a previous step).
    %
    % ---- INUPUT ------------------------------------------------------------
    %       CtrlPts  Control points for Bezier curve, [2x4]
    %         Tvals  Values between 0 an 1, inclusive [1x?]
    %    BezierVals  Points in the Bezier curve [2x?]
    %     CirRadius  Location of the wheel center at values Tvals [2x?]
    %
    % ---- OUTPUT ------------------------------------------------------------
    %   TangentVals  Radius of spirograph wheel, a negative radius indicates
    %                that the wheel rolls inside the curve [1]
    %
    function [BezierNormal] = EvalNormal( obj, TVals, CirRadius)
      % the tangent vector is the derivative of the curve
      BezierPrime = 3*( (1-TVals).^2) .* (obj.CtrlPts(:,2)-obj.CtrlPts(:,1))+ ...
        6 * ((1-TVals).*(TVals)) .* (obj.CtrlPts(:,3)-obj.CtrlPts(:,2)) + ...
        3 * ( (TVals).^2) .* (obj.CtrlPts(:,4)-obj.CtrlPts(:,3)) ;
      % the normal tangent vector is obtained by rotating the tangent vector
      BezierNormal = CirRadius * [0,-1; 1,0]*BezierPrime ./ vecnorm(BezierPrime,2,1);
      
      % interpolation for vanishing derivatives
      if any(vecnorm(BezierPrime,2,1) < 0.001)
        for i = 1:size(TVals,2)
          if norm(BezierPrime(:,i)) < 0.001
            TVals2 = ((-0.001):0.0001:0.001) + TVals(i);
            TVals2 = TVals2( (TVals2~=TVals(i)) & (TVals2>0) & (TVals2<1) );
            BezierPrime2 = 3*( (1-TVals2).^2) .* (obj.CtrlPts(:,2)-obj.CtrlPts(:,1))+ ...
              6 * ((1-TVals2).*(TVals2)) .* (obj.CtrlPts(:,3)-obj.CtrlPts(:,2)) + ...
              3 * ( (TVals2).^2) .* (obj.CtrlPts(:,4)-obj.CtrlPts(:,3)) ;
            BezierNormal2 = CirRadius * [0,-1; 1,0]*BezierPrime2 ./ vecnorm(BezierPrime2,2,1);
            weights = exp( -( ( TVals2-TVals(i) ).^2)/(2*(0.0005)^2) );
            BezierNormal_replacement = BezierNormal2 * weights' / sum(weights);
            BezierNormal(:,i) = BezierNormal_replacement / norm(BezierNormal_replacement);
          end
        end
      end
    end
    % Compute the tangent vector to a cubic Bezier curve at a set of values
    % between 0 and 1.
    % In order to avoid redundant computations, the points on the Bezier curve
    % are taken as arguments (computed on a previous step).
    %
    % ---- INUPUT ------------------------------------------------------------
    %       CtrlPts  Control points for Bezier curve, [2x4]
    %         Tvals  Values between 0 an 1, inclusive [1x?]
    %    BezierVals  Points in the Bezier curve [2x?]
    %     CirRadius  Location of the wheel center at values Tvals [2x?]
    %
    % ---- OUTPUT ------------------------------------------------------------
    %   TangentVals  Radius of spirograph wheel, a negative radius indicates
    %                that the wheel rolls inside the curve [1]
    %
    function [BezierTangent] = EvalBezierTangent( obj, TVals, CirRadius)
      % the tangent vector is the derivative of the curve
      BezierPrime = 3*( (1-TVals).^2) .* (obj.CtrlPts(:,2)-obj.CtrlPts(:,1))+ ...
        6 * ((1-TVals).*(TVals)) .* (obj.CtrlPts(:,3)-obj.CtrlPts(:,2)) + ...
        3 * ( (TVals).^2) .* (obj.CtrlPts(:,4)-obj.CtrlPts(:,3)) ;
      % the normal tangent vector is obtained by rotating the tangent vector
      BezierTangent = CirRadius * BezierPrime ./ vecnorm(BezierPrime,2,1);
      
      % interpolation for vanishing derivatives
      if any(vecnorm(BezierPrime,2,1) < obj.Tol)
        for i = 1:size(TVals,2)
          if norm(BezierPrime(:,i)) < obj.Tol
            TVals2 = ((-10):10)*obj.Tol + TVals(i);
            TVals2 = TVals2( (TVals2~=TVals(i)) & (TVals2>0) & (TVals2<1) );
            BezierPrime2 = 3*( (1-TVals2).^2) .* (obj.CtrlPts(:,2)-obj.CtrlPts(:,1))+ ...
              6 * ((1-TVals2).*(TVals2)) .* (obj.CtrlPts(:,3)-obj.CtrlPts(:,2)) + ...
              3 * ( (TVals2).^2) .* (obj.CtrlPts(:,4)-obj.CtrlPts(:,3)) ;
            BezierNormal2 = CirRadius * [0,-1; 1,0]*BezierPrime2 ./ vecnorm(BezierPrime2,2,1);
            weights = exp( -( ( TVals2-TVals(i) ).^2)/(2*(0.0005)^2) );
            BezierNormal_replacement = BezierNormal2 * weights' / sum(weights);
            BezierTangent(:,i) = BezierNormal_replacement / norm(BezierNormal_replacement);
          end
        end
      end
    end
    %%
    % Corner case: the wheel can't continue rolling on the current curve
    % because it will collide with the next curve.
    % Under such assumption, this function will find the indexes t1, t2 from
    % the curves C1 and C2 such that the scaled normals at C1(t1) and C2(t2) 
    % are equal --up to some tolerance. This means that, it the curve was a
    % physical object, the wheel will stop there.
    % Needless to say, this process depends heavily on the wheel size.
    %
    % ---- INUPUT ------------------------------------------------------------
    %      Curve1/2  Array with control points for each one of the Bezier
    %                curves that make the curve {?} <- [2,4]'s
    %           Tol  Max distance between points of the discretization [1]
    %
    % ---- OUTPUT ------------------------------------------------------------
    %          Time  Timestamps [1x?]
    %      WhCtrPos  Location of the wheel center at timestamps [2x?]
    %     MarkerPos  Location of marker at timepoints [2x?]
    %
    % The last control point of the last curve must be equal to the first
    % control point of the first curve. This is not checked.
    %
    function [CrossTime1, CrossTime2] = ...
      FindCollisionTime( obj, obj2, WheelRadius )
    % initial and finishing times
    t0_1 = 0;
    t0_2 = 0;
    tF_1 = 1;
    tF_2 = 1; 

    % approximate the curve: if precision up to Tol is not reached, iterate
    iter = 1;
    while iter < 50
      % interpolate curves adn their normals
      LocalTime1 = t0_1:((tF_1-t0_1)/(ceil((tF_1-t0_1)/obj.Tol)+iter)):tF_1;
      LocalTime2 = t0_2:((tF_2-t0_2)/(ceil((tF_2-t0_2)/obj.Tol)+iter)):tF_2; 
      %
      BezierPos1  = obj.EvalPosition(  LocalTime1 );
      BezierPos2  = obj2.EvalPosition( LocalTime2 );
      BezierNorm1 = obj.EvalNormal(  LocalTime1, WheelRadius );
      BezierNorm2 = obj2.EvalNormal( LocalTime2, WheelRadius );
      %
      WhCtrPos1 = BezierPos1 + BezierNorm1;
      WhCtrPos2 = BezierPos2 + BezierNorm2;
      %
      % look for the intersection per se
      curve1_dis = zeros(size(LocalTime1));
      curve1_dis_idx = zeros(size(LocalTime1));
        for i = 1:size(LocalTime1,2)
          [dist, idxx] = min( vecnorm( WhCtrPos1(:,i) - WhCtrPos2, 2, 1 ) );
          curve1_dis(i) = dist;
          curve1_dis_idx(i) = idxx;
        end
        %
        % get the intersection time
        [dist, idx1] = min(curve1_dis);
        idx2 = curve1_dis_idx(idx1);
        %
        % if the normals are different, make a better approximation
        if dist < obj.Tol
          break
        end
        t0_1 = max(0,LocalTime1(idx1)-0.5/2^iter);
        t0_2 = max(0,LocalTime2(idx2)-0.5/2^iter);
        tF_1 = min(1,LocalTime1(idx1)+0.5/2^iter);
        tF_2 = min(1,LocalTime2(idx2)+0.5/2^iter);
        iter = iter+1;
      end
      %
      % report results
      if iter < MaxIter
        CrossTime1 = LocalTime1(idx1);
        CrossTime2 = LocalTime2(idx2);
      else
        CrossTime1 = nan;
        CrossTime2 = nan;
      end
    end
  end
end




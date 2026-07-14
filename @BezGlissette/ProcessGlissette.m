function ProcessGlissette( obj )

% this evaluation is for background decoration only
obj.DecorativeBez = obj.BPath.EvalAllPositions();

% additional, specific, preparation for each method
switch obj.Method
  case 'Hole'
    obj.DecorativeHole = obj.HPath.EvalAllPositions();
end

% trying to generate a glissette without defining a method
if ~isprop(obj, 'Method') || isempty(obj.Method)
  obj.Method = 'Default';
end

% default parameters
switch obj.Method
  case 'Default'
    if ~isprop(obj, 'MarkerRadius') || isempty(obj.MarkerRadius)
      warning('Using default values.')
      obj.MarkerRadius = 0;
    end
  case 'Hole'
    if ~isprop(obj, 'BezBase') || ~isprop(obj, 'AngBase')
      warning('Run "SetupHole" before generating glissette.')
      return
    end
  case 'Ring2'
    if ~isprop(obj, 'CtrHoleDist') || ~isprop(obj, 'HoleRadius') || ~isprop(obj, 'Wheel2Radius') || ~isprop(obj, 'Marker2Radius')
      warning('Multiple arguments missing.')
      return
    end
end

% containers for results
Time        = [];
BezierPos   = [];
WhCtrPos    = [];
WhCtrAngle  = [];
MarkerPos   = [];
MarkerAngle = [];

% initialize
CurrAngle0   = obj.MarkerAngle0; % point to the curve
CurrTime0    = 0;

%%
% loop
CurrSpin = 0; % index start at 0
ClosedFlag = false;
SufficientSpins = false;
CurrRollDist0 = 0;
while (CurrSpin < obj.MaxSpins) && (~ClosedFlag)
  for j = 1:obj.BPath.nSegments
    disp(strcat('Spin: ',num2str(CurrSpin),' , Segment: ',num2str(j)))
    CurrSegment = obj.BPath.Segment{j};
    
    % run one single Bezier curve at the time
    % initial guess for time
    LocalTime = linspace(0,1, ceil(CurrSegment.GetSegmentPerimeter()/obj.Tol) );

    % main loop
    iter = 0;
    %locTolFlag = false;
    while (iter < obj.MaxIter) %&& ~locTolFlag
      iter = iter+1; % max iterations
      nPts = size(LocalTime,2);
      
      % compute the points on the Bezier curve and the wheel that rolls over it
      locBezierPos  = CurrSegment.EvalPosition( LocalTime );
      locBezierNorm = CurrSegment.EvalNormal(   LocalTime, obj.Wheel1Radius );
      if obj.Wheel1Radius > 0
        locWh1CtrPos = locBezierPos + locBezierNorm;
      else
        locWh1CtrPos = locBezierPos;
      end
      BezNormAngle = -atan2(locBezierNorm(2,:), locBezierNorm(1,:));

      % distance rolled by wheel 1
      CumRollWh1 = cumsum([CurrRollDist0, vecnorm( diff(locBezierPos,1,2), 2, 1)]);
      if obj.Wheel1Radius > 0
        CumAngleWh1 = CumRollWh1 / obj.Wheel1Radius;
      else
        CumAngleWh1 = CumRollWh1;
      end
      locWh1Angle = -CumAngleWh1 - BezNormAngle;

      % position and angle for the marker
      switch obj.Method
        case 'Default'
          % add a point in a circle with given center and radius
          locMarkerPos = locWh1CtrPos + [cos(locWh1Angle); sin(locWh1Angle)]*obj.MarkerRadius;
        case 'Hole'
          % distance that the wheel has rolled so far
          % 1. interpolate where the marker is in the hole
          HolePos = zeros(2, nPts);
          HolePos(1,:) = interp1(obj.AngBase, obj.BezBase(1,:), mod(CumAngleWh1,2*pi));
          HolePos(2,:) = interp1(obj.AngBase, obj.BezBase(2,:), mod(CumAngleWh1,2*pi));
          % 2. rotate the interpolated point and add around the wheel center
          locMarkerPos = zeros(2, nPts);
          for k = 1:nPts
            th = locWh1Angle(k)+pi;
            locMarkerPos(:,k) = locWh1CtrPos(:,k) + ...
              [cos(th) -sin(th); sin(th) cos(th)] * HolePos(:,k) * obj.Wheel1Radius;
          end
        case 'Ring2'
          % distance that the wheel has rolled so far
          % compute the angle rolled by the smallest gear
          CircProject = [cos(CumAngleWh1); sin(CumAngleWh1)] - [obj.CtrHoleDist;0];
          LocRollAngle = atan2(CircProject(2,:),CircProject(1,:));
          % position of marker IF ring was static
          InnMarkerPos = [...
            (HoleRadius-Wheel2Radius)*cos(LocRollAngle) + Marker2Radius*cos(LocRollAngle*((HoleRadius-Wheel2Radius)/Wheel2Radius));...
            (HoleRadius-Wheel2Radius)*sin(LocRollAngle) - Marker2Radius*sin(LocRollAngle*((HoleRadius-Wheel2Radius)/Wheel2Radius))...
            ];
          % rotate the interpolated point and add around the wheel center
          locMarkerPos = zeros(2, nPts);
          for k = 1:nPts
            th = locWh1Angle(k);
            locMarkerPos(:,k) = locWh1CtrPos(:,k) + ...
              [cos(th) -sin(th); sin(th) cos(th)] * ( InnMarkerPos(:,k) + [obj.CtrHoleDist;0] );
          end
        otherwise
          locMarkerPos = locWh1CtrPos;
      end
      
      % check if the marker points are not too far from each other
      DiffCurve = vecnorm( diff(locMarkerPos,1,2), 2, 1);
      if max(DiffCurve) < obj.Tol
        %locTolFlag = true;
        break
      end
      
      % if the marker points are too far, add more time points when needed
      NewTimes = [];
      for i = 2:length(LocalTime)
        if( DiffCurve(i-1) > obj.Tol )
          epsilon = (LocalTime(i)-LocalTime(i-1))/ceil(DiffCurve(i-1)/(obj.Tol/2));
          NewTimes = [ NewTimes, (LocalTime(i-1):epsilon:LocalTime(i)) ];
        end
      end

      % use the new timepoints and iterate until the result is acceptable
      LocalTime = unique( [LocalTime, NewTimes], "sorted" );
    end

    % technical
    locWh1Angle = mod(locWh1Angle, 2*pi);

    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        LocalTime + CurrTime0];
    BezierPos   = [BezierPos,   locBezierPos  ];
    WhCtrPos    = [WhCtrPos,    locWh1CtrPos  ];
    WhCtrAngle  = [WhCtrAngle,  mod(CumAngleWh1,2*pi)];
    MarkerPos   = [MarkerPos,   locMarkerPos  ];
    MarkerAngle = [MarkerAngle, locWh1Angle   ];

    %
    % update initial values
    CurrTime0  = Time(end);
    CurrRollDist0 = CurrRollDist0 + CurrSegment.GetSegmentPerimeter();

    % debug
    if false
      figure()
      axis equal
      xlim([-1,1])
      ylim([-1,1])
      hold on
      grid on
      plot(HolePos(1,:), HolePos(2,:))
      scatter(HolePos(1,end),HolePos(2,end),'filled')
    end

    % debug
    if false
      figure()
      hold on
      axis equal
      grid on
      xlim([-1,1])
      ylim([-1,1])
      plot(BezierPos(1,:),BezierPos(2,:))
      scatter(BezierPos(1,end),BezierPos(2,end),'filled')
      plot(WhCtrPos(1,:),WhCtrPos(2,:))
      scatter(WhCtrPos(1,end),WhCtrPos(2,end),'filled')
      plot(MarkerPos(1,:),MarkerPos(2,:))
      scatter(MarkerPos(1,end),MarkerPos(2,end),'filled')
    end
    
    % prepare for a corner
    NextSegment = obj.BPath.Segment{mod(j+1-1,obj.BPath.nSegments)+1};
    
    % roll over the corner, if needed
    %%
    % compute arc angle that the wheel center will describe
    WhNormalPre = CurrSegment.EvalNormal(1, obj.Wheel1Radius);
    WhNormalPos = NextSegment.EvalNormal(0, obj.Wheel1Radius);

    % early stop if the wheel won't actually roll
    if (norm( WhNormalPos-WhNormalPre ) < obj.Tol)&&(norm( NextSegment.EvalPosition(0)-CurrSegment.EvalPosition(1) ) < obj.Tol)
      continue
    end

    % early stop is found a discontinuity
    if (norm( NextSegment.EvalPosition(0)-CurrSegment.EvalPosition(1) ) < obj.Tol)
      % TO BE PATCHED
      continue
    end

    CornerAngle = atan2(WhNormalPos(2),WhNormalPos(1)) - atan2(WhNormalPre(2),WhNormalPre(1));
    if CornerAngle > 0
      CornerAngle = -(2*pi-CornerAngle);
    end

    % compute the length that the marker will describe
    MarkerPos0 = MarkerPos(:,end);
    BezierPos0 = BezierPos(:,end);
    MarkerArcLength = abs(CornerAngle) * norm( MarkerPos0 - BezierPos0 );

    % adjust the rotation, time is set in terms of the marker
    LocalTime   = linspace(0,1, ceil(MarkerArcLength / obj.Tol));
    locWh1Angle = MarkerAngle0 + LocalTime*CornerAngle;

    LocalWhCtrAngle   = atan2(WhNormalPre(2),WhNormalPre(1));
    LocalMarkerAngle  = atan2(MarkerPos0(2)-BezierPos0(2),MarkerPos0(1)-BezierPos0(1));
    LocalMarkerRadius = norm( MarkerPos0 - BezierPos0 );

    locWh1CtrPos  = BezierPos0 + [cos(LocalWhCtrAngle  + LocalTime*CornerAngle); sin(LocalWhCtrAngle  + LocalTime*CornerAngle)]*obj.Wheel1Radius;
    locMarkerPos = BezierPos0 + [cos(LocalMarkerAngle + LocalTime*CornerAngle); sin(LocalMarkerAngle + LocalTime*CornerAngle)]*LocalMarkerRadius;

    % technical
    locWh1Angle = mod(locWh1Angle, 2*pi);
    BezierPos   = BezierPos0 * ones(size(LocalTime));

    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        LocalTime+CurrTime0];
    BezierPos   = [BezierPos,   locBezierPos];
    WhCtrPos    = [WhCtrPos,    locWh1CtrPos];
    MarkerPos   = [MarkerPos,   locMarkerPos];
    MarkerAngle = [MarkerAngle, locWh1Angle];

    % update initial values
    CurrTime0  = Time(end);
  end
  %
  %%
  % update number of spins
  CurrSpin = CurrSpin + 1;
  if CurrSpin >= obj.MinSpins
    SufficientSpins = true;
  end
  %
  % check if spirograph is closed
  if CurrSpin == 1
    FirstPt = MarkerPos(:,1);
  end
  LastPt = MarkerPos(:,end);
  if norm( FirstPt - LastPt ) < obj.CloseTol
    if SufficientSpins
      ClosedFlag = true;
    end
  end
end

%%

% patch
if obj.CloseEnds
  MarkerPos(:,end+1) = MarkerPos(:,1);
  MarkerAngle(end+1) = MarkerAngle(1);
end

% reverse curves with opposite orientation
if isprop(obj,'ChangeOrient') && obj.ChangeOrient
  BezierPos   = flip(BezierPos, 2);
  WhCtrPos    = flip( WhCtrPos, 2 );
  MarkerPos   = flip(MarkerPos, 2);
  MarkerAngle = flip(MarkerAngle, 2);
  Time = max(Time) - flip(Time);
end

% DEBUG: only assign after all process are complete
obj.BezierPos   = BezierPos;
obj.WhCtrPos    = WhCtrPos;
obj.MarkerPos   = MarkerPos;
obj.MarkerAngle = MarkerAngle;
obj.LocTime     = Time;

end
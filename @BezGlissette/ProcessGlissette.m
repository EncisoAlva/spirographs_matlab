function ProcessGlissette( obj )

% this evaluation is for background decoration only
obj.DecorativeBez = obj.BPath.EvalPotition( obj.Tol );

% trying to generate a glissette without defining a method
if ~isprop(obj, 'Method')
  obj.Method = 'Default';
end

% default parameters
switch obj.Method
  case 'Default'
    if ~isprop(obj, 'MarkerRadius')
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
MarkerPos   = [];
MarkerAngle = [];

% initialize
FirstTangent = obj.BPath.Segment{1}.EvalNormal(0,1);
CurrAngle0   = atan2(FirstTangent(2), FirstTangent(1)) + MarkerAngle0 + pi; % point to the curve
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
    CurrSegment = obj.BPath{j};
    CurrCtrlPts = CurrSegment.CtrlPts;
    
    % run one single Bezier curve at the time
    % initial guess for time
    PerUpBound = ...
      norm(CurrCtrlPts(:,1)-CurrCtrlPts(:,2)) + ...
      norm(CurrCtrlPts(:,2)-CurrCtrlPts(:,3)) + ...
      norm(CurrCtrlPts(:,3)-CurrCtrlPts(:,4));
    DistDelta = 1/ceil(PerUpBound/Tol);
    LocalTime = 0:DistDelta:1;

    % main loop
    iter = 0;
    locTolFlag = false;
    while (iter < obj.MaxIter) && ~locTolFlag
      iter = iter+1; % max iterations
      nPts = size(LocalTime,2);
      
      % compute the points on the Bezier curve and the wheel that rolls over it
      locBezierPos  = CurrSegment.EvalPosition( LocalTime );
      locBezierNorm = CurrSegment.EvalNormal(   LocalTime, obj.Wheel1Radius );
      if obj.Wheel1Radius > 0
        locWhCtrPos = locBezierPos + locBezierNorm;
      else
        locWhCtrPos = locBezierPos;
      end
      BezNormAngleDiff = diff( atan2(locBezierNorm(2,:), locBezierNorm(1,:)), 1, 2 );

      % angle that the wheel spun between two given points
      if obj.Wheel1Radius > 0
        DiffAngle = vecnorm( diff(locBezierPos,1,2), 2, 1) / obj.Wheel1Radius;
      else
        DiffAngle = vecnorm( diff(locBezierPos,1,2), 2, 1);
      end

      % position and angle for the marker
      locMarkerAngle = cumsum([CurrAngle0, -DiffAngle+BezNormAngleDiff]);
      switch obj.Method
        case 'Default'
          % add a point in a circle with given center and radius
          locMarkerPos = locWhCtrPos + [cos(locMarkerAngle); sin(locMarkerAngle)]*obj.MarkerRadius;
        case 'Hole'
          % distance that the wheel has rolled so far
          RollAngle = mod( cumsum([RollDist0, DiffAngle ]), 2*pi);
          % 1. interpolate where the marker is in the hole
          HolePos = zeros(2, nPts);
          HolePos(1,:) = interp1(obj.AngBase, obj.BezBase(1,:), mod(RollAngle,2*pi));
          HolePos(2,:) = interp1(obj.AngBase, obj.BezBase(2,:), mod(RollAngle,2*pi));
          % 2. rotate the interpolated point and add around the wheel center
          locMarkerPos = zeros(2, nPts);
          for k = 1:nPts
            th = locMarkerAngle(k);
            locMarkerPos(:,k) = locWhCtrPos(:,k) + ...
              obj.Wheel1Radius * [cos(th) -sin(th); sin(th) cos(th)] * HolePos(:,j);
          end
        case 'Ring2'
          % distance that the wheel has rolled so far
          RollAngle = cumsum([RollDist0, DiffAngle ]);
          % compute the angle rolled by the smallest gear
          CircProject = [cos(RollAngle); sin(RollAngle)] - [obj.CtrHoleDist;0];
          LocRollAngle = atan2(CircProject(2,:),CircProject(1,:));
          % position of marker IF ring was static
          InnMarkerPos = [...
            (HoleRadius-Wheel2Radius)*cos(LocRollAngle) + Marker2Radius*cos(LocRollAngle*((HoleRadius-Wheel2Radius)/Wheel2Radius));...
            (HoleRadius-Wheel2Radius)*sin(LocRollAngle) - Marker2Radius*sin(LocRollAngle*((HoleRadius-Wheel2Radius)/Wheel2Radius))...
            ];
          % rotate the interpolated point and add around the wheel center
          locMarkerPos = zeros(2, nPts);
          for k = 1:nPts
            th = locMarkerAngle(k);
            locMarkerPos(:,k) = locWhCtrPos(:,k) + ...
              [cos(th) -sin(th); sin(th) cos(th)] * ( InnMarkerPos(:,k) + [obj.CtrHoleDist;0] );
          end
        otherwise
          locMarkerPos = locWhCtrPos;
      end
      
      % check if the marker points are not too far from each other
      DiffCurve = vecnorm( diff(MarkerPos,1,2), 2, 1);
      if max(DiffCurve) < Tol
        locTolFlag = true;
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
    locMarkerAngle = mod(locMarkerAngle, 2*pi);

    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        locTime+CurrTime0];
    BezierPos   = [BezierPos,   locBezierPos];
    WhCtrPos    = [WhCtrPos,    locWhCtrPos];
    MarkerPos   = [MarkerPos,   locMarkerPos];
    MarkerAngle = [MarkerAngle, locMarkerAngle];
    %
    % update initial values
    CurrTime0  = Time(end);
    CurrAngle0 = MarkerAngle(end);
    CurrRollDist0 = CurrRollDist0 + PathPerimeter(CurrCtrlPts, Tol)/obj.Wheel1Radius;
    
    % prepare for a corner
    NextSegment = obj.BPath{mod(j+1-1,obj.BPath.nSegments)+1};
    NextCtrlPts = NextSegment.CtrlPts;
    
    % roll over the corner, if needed
    %%
    % compute arc angle that the wheel center will describe
    WhNormalPre = CurrSegment.EvalNormal(1, obj.Wheel1Radius);
    WhNormalPos = NextSegment.EvalNormal(0, obj.Wheel1Radius);

    % early stop if the wheel won't actually roll
    if norm( (NextCtrlPts(:,1)+WhNormalPos)-(CurrCtrlPts(:,4)+WhNormalPre) ) < obj.Tol
      locTime = [];
      locBezierPos = [];
      locWhCtrPos = [];
      locMarkerPos = [];
      locMarkerAngle = [];
      continue
    end

    CornerAngle = atan2(WhNormalPos(2),WhNormalPos(1)) - atan2(WhNormalPre(2),WhNormalPre(1));
    if CornerAngle > 0
      CornerAngle = -(2*pi-CornerAngle);
    end

    % compute the length that the marker will describe
    MarkerPos0 = MarkerPos(:,end);
    MarkerArcLength = abs(CornerAngle) * norm( MarkerPos0 - NextCtrlPts(:,1) );

    % adjust the rotation, time is set in terms of the marker
    LocalTime   = 0:1/ceil(MarkerArcLength / obj.Tol):1;
    locMarkerAngle = MarkerAngle0 + LocalTime*CornerAngle;

    LocalWhCtrAngle   = atan2(WhNormalPre(2),WhNormalPre(1));
    LocalMarkerAngle  = atan2(MarkerPos0(2)-CurrCtrlPts(2,end),MarkerPos0(1)-CurrCtrlPts(1,end));
    LocalMarkerRadius = norm( MarkerPos0 - NextCtrlPts(:,1) );

    locWhCtrPos  = CurrCtrlPts(:,4) + [cos(LocalWhCtrAngle  + LocalTime*CornerAngle); sin(LocalWhCtrAngle  + LocalTime*CornerAngle)]*obj.Wheel1Radius;
    locMarkerPos = CurrCtrlPts(:,4) + [cos(LocalMarkerAngle + LocalTime*CornerAngle); sin(LocalMarkerAngle + LocalTime*CornerAngle)]*LocalMarkerRadius;

    locTime = LocalTime + Time0;

    % technical
    locMarkerAngle = mod(locMarkerAngle, 2*pi);
    BezierPos   = CurrCtrlPts(:,1) * ones(size(locTime));

    %
    % concatenate results from the current segment to the overall outputs
    Time        = [Time,        locTime+CurrTime0];
    BezierPos   = [BezierPos,   locBezierPos];
    WhCtrPos    = [WhCtrPos,    locWhCtrPos];
    MarkerPos   = [MarkerPos,   locMarkerPos];
    MarkerAngle = [MarkerAngle, locMarkerAngle];

    % update initial values
    CurrTime0  = Time(end);
    CurrAngle0 = MarkerAngle(end);
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
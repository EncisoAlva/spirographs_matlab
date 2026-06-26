% Take a beier path inside an unit circle and simulate the rotation of the
% circle as a point travels through the path. 
% This will result on a table of angle vs position-at-path, which willbe
% interpolated later to simulate the rotation of said circle over another
% curve.
%
% The relationship between path transversed on the path vs angle rotated is
% ASSUMED to be monotopnic; a sufficient condition is that the bezier path 
% is convex.
%
% ---- INUPUT ------------------------------------------------------------
%         HPath  Bezier path representing the hole in the rolling wheel.
%           Tol  Maximum allowable distance between neighboring points [1]
%      CheckFit  (Bool) Show a graph with the resulting associations [1]
%
% ---- OUTPUT ------------------------------------------------------------
%       BezBase  Points B(p) at the bezier path [2x?]
%       AngBase  Angles A(p) tied to B(p) [1x?]
%       DisBase  Distance p(t) = lineint_0^t B(t) dt
%
function SetupHole( obj, CheckFit )

% find values of t with a prescribed max dist between them
allBez = cell(1,obj.nSegments*2);
allAng = cell(1,obj.nSegments*2);

for s = 1:obj.nSegments
  currSegment = obj.HPath.Segment{s};
  % initial guess for tvals
  Tvals = linspace(0,1, ceil(currSegment.GetSegmentPerimeter()/obj.Tol)+1);
  nTPts = size(Tvals,2);
  %
  tol_flag = false;
  while ~tol_flag
    Bez =  currSegment.EvalPosition(Tvals);
    Nor = -currSegment.EvalNormal(  Tvals, 1);
    %
    % current side
    % reach circumference
    Ksc = zeros(1, nTPts);
    BezProj = zeros(1, nTPts);
    for i = 1:nTPts
      Ksc(i)  = -Bez(:,i)'*Nor(:,i) + sqrt((Bez(:,i)'*Nor(:,i))^2 +1 -norm(Bez(:,i))^2);
      BezProj(:,i) = Bez(:,i) + Ksc(i)*Nor(:,i);
    end
    Ang = atan2(BezProj(2,:),BezProj(1,:));
    %
    % error correction to be added
    tol_flag = true;
  end
  %
  % report results
  allBez{2*s-1} = Bez;
  allAng{2*s-1} = Ang;
  %
  % between sides, only if necessary
  if s<obj.nSegments
    nextSegment = obj.HPath.Segment{s+1};
  else
    nextSegment = obj.HPath.Segment{1};
  end
  currNor = -currSegment.EvalNormal( 1, 1);
  nextNor = -nextSegment.EvalNormal( 0, 1);
  if atan2(nextNor(2),nextNor(1)) - atan2(currNor(2),currNor(1)) ~= 0
  Bez = currSegment.EvalPosition(1);
  if abs(norm(Bez)-1) > obj.Tol % if not currently touching the circle
    % doing this exactly twice, so I am recycling variables
    %Bez = EvalBezier(HPath{s},1);
    Ksc = -Bez'*currNor + sqrt((Bez'*currNor)^2 +1 -norm(Bez)^2);
    BezProj = Bez + Ksc*currNor;
    currAng = mod( atan2(BezProj(2),BezProj(1)), 2*pi);
    %
    Bez = nextSegment.EvalPosition(0);
    Ksc = -Bez'*nextNor + sqrt((Bez'*nextNor)^2 +1 -norm(Bez)^2);
    BezProj = Bez + Ksc*nextNor;
    nextAng = mod( atan2(BezProj(2),BezProj(1)), 2*pi);
    %
    Ang = linspace(currAng, nextAng, ceil((nextAng-currAng)/obj.Tol)+1);
    %
    % report reults
    allBez{2*s} = Bez*ones(size(Ang)); % repeat as needed
    allAng{2*s} = Ang;
  end
  end
end

% move from cell array to vector
Bezz = [];
Angg = [];
for ii = 1:(2*obj.nSegments)
  Bezz = [Bezz, allBez{ii}];
  Angg = [Angg, allAng{ii}];
end

% sort
Angg = mod(Angg, 2*pi);
[Angg, idx] = sort(Angg);
Bezz = Bezz(:,idx);

% remove duplicates
[~, i, ~] = unique(Angg,'first');
idxDupes  = find(not(ismember(1:numel(Angg),i)));
%
Angg(   idxDupes) = [];
Bezz( :,idxDupes) = [];

% if either 0 or 2*pi are not accounted for
if ~ismember(0, Angg) || ~ismember(2*pi, Angg)
  % pick angles around 0 to interpolate
  idx0 = 1:length(Angg);
  idxA = idx0(abs(Angg     )<0.2);
  idxB = idx0(abs(Angg-2*pi)<0.2);
  idxx = unique([idxA, idxB]);
  if ~isempty(idxx)
  Ang_int = Angg(    idxx );
  Bez_int = Bezz( :, idxx );
  Ang_int = mod(Ang_int +pi, 2*pi) -pi; % find equivalents from ~2pi to ~0
  % interpolate
  Bez0 = zeros(2,1);
  Bez0(1) = interp1(Ang_int, Bez_int(1,:), 0);
  Bez0(2) = interp1(Ang_int, Bez_int(2,:), 0);
  if ~ismember(0, Angg)
    Bezz = [Bezz, Bez0];
    Angg = [Angg, 0];
  end
  if ~ismember(2*pi, Angg)
    Bezz = [Bezz, Bez0];
    Angg = [Angg, 2*pi];
  end
  end
end

% report results
obj.BezBase = Bezz;
obj.AngBase = Angg;

if CheckFit
  figure()
  hold on
  xlim([-1,1])
  ylim([-1,1])
  axis equal
  for jj = 1:size(Bezz,2)
    plot([Bezz(1,jj), cos(Angg(jj))], [Bezz(2,jj), sin(Angg(jj))],'blue')
  end
end

end
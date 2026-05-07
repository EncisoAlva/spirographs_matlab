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
function [BezBase, AngBase, DisBase] = SetupHole(HPath, Tol, CheckFit)

% find values of t with a prescribed max dist between them
nSegments = size(HPath,2);
allBez    = cell(1,nSegments*2);
allAng    = cell(1,nSegments*2);

for s = 1:nSegments
  currSegment = HPath{s};
  % initial guess for tvals
  Tvals = linspace(0,1, ceil(PathPerimeter({currSegment},Tol)/Tol)+1);
  %
  tol_flag = false;
  while ~tol_flag
    Bez = EvalBezier(currSegment, Tvals);
    Nor = -EvalBezierNormal(currSegment, Tvals, 1);
    %
    % current side
    % reach circumference
    Ksc = zeros(size(Tvals));
    BezProj = zeros(size(Bez));
    for i = 1:size(Tvals,2)
      Ksc(i)  = -Bez(:,i)'*Nor(:,i) + sqrt((Bez(:,i)'*Nor(:,i))^2 +1 -norm(Bez(:,i))^2);
      BezProj(:,i) = Bez(:,i) + Ksc(i)*Nor(:,i);
    end
    Ang = atan2(BezProj(2,:),BezProj(1,:));
    %
    % error correction to be added
    tol_flag = true;
  end
  %
  % report reults
  allBez{2*s-1} = Bez;
  allAng{2*s-1} = Ang;
  %
  % between sides, only if necessary
  if s<nSegments
    s_next = s+1;
  else
    s_next = 1;
  end
  currNor = -EvalBezierNormal(HPath{s}, 1, 1);
  nextNor = -EvalBezierNormal(HPath{s_next}, 0, 1);
  if atan2(nextNor(2),nextNor(1)) - atan2(currNor(2),currNor(1)) ~= 0
    % doing this exactly twice, so I am recycling variables
    Bez = EvalBezier(HPath{s},1);
    Ksc = -Bez'*currNor + sqrt((Bez'*currNor)^2 +1 -norm(Bez)^2);
    BezProj = Bez + Ksc*currNor;
    currAng = mod( atan2(BezProj(2),BezProj(1)), 2*pi);
    %
    Bez = EvalBezier(HPath{s_next},0);
    Ksc = -Bez'*nextNor + sqrt((Bez'*nextNor)^2 +1 -norm(Bez)^2);
    BezProj = Bez + Ksc*nextNor;
    nextAng = mod( atan2(BezProj(2),BezProj(1)), 2*pi);
    %
    Ang = linspace(currAng, nextAng, ceil((nextAng-currAng)/Tol)+1);
    %
    % report reults
    allBez{2*s} = Bez*ones(size(Ang)); % repeat as needed
    allAng{2*s} = Ang;
  end
end

% move from cell array to vector
Bezz = [];
Angg = [];
for ii = 1:(2*nSegments)
  Bezz = [Bezz, allBez{ii}];
  Angg = [Angg, allAng{ii}];
end
% for convenience, the rotation will start at 0
Bezz = [Bezz, [1;0]];
Angg = [Angg, 0];

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

% PATCH: adding the duplicate angle 0=2pi for interpolation
Bezz = [Bezz, [1;0]];
Angg = [Angg, 2*pi];

% report results
BezBase = Bezz;
AngBase = Angg;

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

% fix later
DisBase = 0;

end
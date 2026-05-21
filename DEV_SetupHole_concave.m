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
function [BezBase, AngBase, DisBase] = DEV_SetupHole_concave(HPath, Tol, CheckFit)

% find values of t with a prescribed max dist between them
nSegments = size(HPath,2);

% will replace with tables but I do not have the time right now
% containers
tab_Bez  = [];
tab_T    = [];
tab_Segm = [];

% identify convex hull
for s = 1:nSegments
  currSegment = HPath{s};
  % initial guess for tvals
  Tvals = linspace(0,1, ceil(PathPerimeter({currSegment},Tol)/Tol)+1);
  %
  tol_reached = false;
  while ~tol_reached
    Bez = EvalBezier(currSegment, Tvals);
    %
    % add later error correction
    tol_reached = true;
  end
  tab_Bez = [tab_Bez, Bez];
  tab_T   = [tab_T, Tvals];
  tab_Segm = [tab_Segm, s*ones(1,size(Tvals,2))];
end

% index of points in convex hull
CHind = convhull(tab_Bez')';

% add duplicates on the convex hull
for j = 2:size(tab_T,2)
  if ismember(j,CHind)
    if norm(tab_Bez(:,j)-tab_Bez(:,j-1)) == 0
      % if a repeated point
      CHind = [CHind, j];
    end
  end
end
if ismember(1,CHind)
  if norm(tab_Bez(:,end)-tab_Bez(:,1)) == 0
    CHind = [CHind, size(tab_T,2)];
  end
end
CHind = sort(unique(CHind));

% identify 'gaps' in the hole
GapIdx = [];
for i = 2:size(CHind,2)
  if CHind(i)-CHind(i-1) > 1
    GapIdx = [GapIdx, [CHind(i-1); CHind(i)]];
  end
end
nGaps = size(GapIdx,2);

% label points according to the gap
tab_Gap = zeros(1,size(tab_T,2));
currGap = 1;
ingap = false;
for i = 1:size(tab_T,2)
  if ~ingap
    if ismember(i,CHind)
      % no change, no gap -> no gap
    else
      % entering gap
      ingap = true;
      tab_Gap(i) = currGap;
    end
  else
    if ismember(i,CHind)
      % getting outside the gap
      ingap = false;
      currGap = currGap + 1;
    else
      % no change, gap -> gap
      tab_Gap(i) = currGap;
    end
  end
end

% line points in gaps
GapStart = zeros(2,nGaps);
GapVect  = zeros(2,nGaps);
for j = 1:nGaps
  GapStart(:,j) = tab_Bez(:,GapIdx(1,j));
  tmp = tab_Bez(:,GapIdx(2,j)) - tab_Bez(:,GapIdx(1,j));
  GapVect(:,j) = tmp/norm(tmp);
end

% compute tangents
Ang = zeros(1,size(tab_T,2));
for i = 1:size(tab_T,2)
  if ismember(i,CHind)
    % normal vector
    Bez_i = tab_Bez(:,i);
    Nor_i = -EvalBezierNormal(HPath{tab_Segm(i)}, tab_T(i), 1);
    % reach circumference
    Ksc  = -Bez_i'*Nor_i + sqrt((Bez_i'*Nor_i)^2 +1 -norm(Bez_i)^2);
    BezProj = Bez_i + Ksc*Nor_i;
    Ang(i) = atan2(BezProj(2),BezProj(1));
  else
    Bez_i = tab_Bez(:,i);
    LinePt = GapStart(:,tab_Gap(i)) + (GapVect(:,tab_Gap(i))'*(Bez_i-GapStart(:,tab_Gap(i))))*GapVect(:,tab_Gap(i));
    tmp = LinePt - Bez_i;
    LineVect = tmp / norm(tmp);
    % reach circumference
    Ksc  = -LinePt'*LineVect + sqrt((LinePt'*LineVect)^2 +1 -norm(LinePt)^2);
    BezProj = LinePt + Ksc*LineVect;
    Ang(i) = atan2(BezProj(2),BezProj(1));
  end
end

% patch
Angg = Ang;
Bezz = tab_Bez;

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
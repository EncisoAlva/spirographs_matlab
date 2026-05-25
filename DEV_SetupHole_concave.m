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
CHind = sort(unique(CHind));

% for any segment, each endpoint is also in another segment
% add duplicates on the convex hull, which is avoided by default
for j = 1:size(tab_T,2)
  if ismember(j,CHind)
    j_ = mod(j-1-1,size(tab_T,2))+1;
    if norm(tab_Bez(:,j)-tab_Bez(:,j_)) == 0
      j_m = mod(j-2-1,size(tab_T,2))+1;
      j_p = mod(j+1-1,size(tab_T,2))+1;
      if ~ismember(j_m,CHind) || ~ismember(j_p,CHind)
        % at border of gap, consider carefully
        if ~ismember(j_m,CHind)
          % gap -> these points -> CH
          CHind(CHind==j_) = [];
          CHind = [CHind, j];
        else
          % CH -> these points -> gap
          CHind(CHind==j) = [];
          CHind = [CHind, j_];
        end
      else
        % away from gaps, just add it
        CHind = [CHind, j_];
      end
    end
    j_ = mod(j+1-1,size(tab_T,2))+1;
    if norm(tab_Bez(:,j)-tab_Bez(:,j_)) == 0
      j_m = mod(j-1-1,size(tab_T,2))+1;
      j_p = mod(j+2-1,size(tab_T,2))+1;
      if ~ismember(j_m,CHind) || ~ismember(j_p,CHind)
        % at border of gap, consider carefully
        if ~ismember(j_m,CHind)
          % gap -> these points -> CH
          CHind(CHind==j) = [];
          CHind = [CHind, j_];
        else
          % CH -> these points -> gap
          CHind(CHind==j_) = [];
          CHind = [CHind, j];
        end
      else
        % away from gaps, just add it
        CHind = [CHind, j_];
      end
    end
  end
end
CHind = sort(unique(CHind));

% identify 'gaps' in the hole
GapIdx = [];
for i = 1:size(CHind,2)
  i_ = mod(i+1-1,size(CHind,2)) +1;
  if abs(CHind(i)-CHind(i_)) > 1
    GapIdx = [GapIdx, [CHind(i); CHind(i_)]];
  end
end
nGaps = size(GapIdx,2);

% label points according to the gap
tab_Gap = zeros(1,size(tab_T,2));
for currGap = 1:nGaps
  i = mod(GapIdx(1,currGap)+1-1,size(tab_T,2))+1;
  while i ~= GapIdx(2,currGap)
    tab_Gap(i) = currGap;
    i = mod(i+1-1,size(tab_T,2))+1;
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
    %tmp = LinePt - Bez_i;
    %LineVect = tmp / norm(tmp);
    LineVect = [0,1; -1,0] * GapVect(:,tab_Gap(i));
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
% patch: what if there is an angle similar to zero, up to tol?
dist_to_0   = min(abs(Angg-0));
dist_to_2pi = min(abs(Angg-2*pi));
if (dist_to_0>Tol/10) || (dist_to_2pi>Tol/10)
  % if at least one is present, copy-paste to the other
  if (dist_to_0<Tol/10) || (dist_to_2pi<Tol/10)
    if (dist_to_0<Tol/10)
      [~,idx_tmp]   = min(abs(Angg-0));
      Bezz = [Bezz, Bezz(:,idx_tmp)];
      Angg = [Angg, 2*pi];
    else
      [~,idx_tmp]   = min(abs(Angg-2*pi));
      Bezz = [Bezz, Bezz(:,idx_tmp)];
      Angg = [Angg, 0];
    end
  else
    % pick angles around 0 to interpolate
    idx0 = 1:length(Angg);
    idxA = idx0(abs(Angg     )<pi/4);
    idxB = idx0(abs(Angg-2*pi)<pi/4);
    idxx = unique([idxB, idxA]);
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
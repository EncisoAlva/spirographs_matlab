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

% ratio between perimeters
PerimRatio = 2*pi/ PathPerimeter(HPath, 0.01);

% find values of t with a prescribed max dist between them
nSegments = size(HPath,2);
allTvals  = cell(1,nSegments);
allBez    = cell(1,nSegments);
allAng    = cell(1,nSegments);
allNor    = cell(1,nSegments);
allTan    = cell(1,nSegments);

currAng = 0;
for ii = 1:nSegments
  % finding values of t
  currSegment = HPath{ii};
  nn = ceil( sum( vecnorm( diff(currSegment,1,2), 2, 1 ) )/Tol );
  Tvals  = linspace(0,1, nn);
  %
  % error correction to be added
  %
  % compute positions B(t_i) and distances |B(t_i)-B(t_i+1)|
  Bez = EvalBezier(currSegment, Tvals);
  Nor = EvalBezierNormal(currSegment, Tvals, 1);
  Tan = EvalBezierTangent(currSegment, Tvals, 1);
  Dis = vecnorm( diff(Bez, 1, 2), 2, 1);
  %
  % cumulative angle
  Ang = [0, cumsum(Dis)]*PerimRatio + currAng;
  currAng = Ang(end);
  %
  % combine results
  allTvals{ii} = Tvals;
  allBez{ii}   = Bez;
  allAng{ii}   = Ang;
  allNor{ii}   = Nor;
  allTan{ii}   = Tan;
end
% patch: ignore rounding errors from adding angles
for ii = 1:nSegments
  allAng{ii} = allAng{ii} * 2*pi/currAng;
end

% move from cell array to vector
Bezz = [];
Angg = [];
Norr = [];
Tann = [];
for ii = 1:nSegments
  Bezz = [Bezz, allBez{ii}];
  Angg = [Angg, allAng{ii}];
  Norr = [Norr, allNor{ii}];
  Tann = [Tann, allTan{ii}];
end

% remove duplicates
[~, i, ~] = unique(Angg,'first');
idxDupes  = find(not(ismember(1:numel(Angg),i)));
%
Angg(   idxDupes) = [];
Bezz( :,idxDupes) = [];
Norr( :,idxDupes) = [];
Tann( :,idxDupes) = [];

% iterative correction: angle increment PROPORTIONAL to area covered
Angg_old = Inf(size(Angg));
while max(abs(Angg - Angg_old)) > 0.01
  %disp(max(abs(Angg - Angg_old)))
  Angg_old = Angg;
  %
  % QUANT ~ area of B(t_i)-C(t_i)-C(t_i+1)-B(t_i+1) 
  DiffAngg     = [0, diff(Angg)];
  DistBezzCirc = vecnorm( Bezz - [cos(Angg);sin(Angg)], 2, 1);
  Discrepancy = zeros(1, size(Bezz,2));
  for jj = 1:size(Bezz,2)
    vecN = (Bezz(:,jj) - [cos(Angg(jj));sin(Angg(jj))]) / norm((Bezz(:,jj) - [cos(Angg(jj));sin(Angg(jj))]));
    Discrepancy(jj) = abs( vecN' * (-Norr(:,jj)) );
  end
  Pull = zeros(1, size(Bezz,2));
  for jj = 1:size(Bezz,2)
    Pull(jj) = abs( [-sin(Angg(jj));cos(Angg(jj))]' * (Tann(:,jj)) );
  end
  %
  %QUANT1 = ( DiffAngg.*DistBezzCirc/2 ).^0.25;
  medianDist = median(DistBezzCirc);
  %QUANT1 = ( DistBezzCirc/medianDist ).^2;
  %QUANT2 = 1-Discrepancy;
  QUANT1 = Pull;
  %QUANT = (((Discrepancy)+0)/1).* DistBezzCirc;
  %Angg  = (Angg + cumsum(QUANT1) * 2*pi/sum(QUANT1) + (cumsum(QUANT2) ) * 2*pi/sum(QUANT2) )/3;
  Angg  = (Angg + cumsum(QUANT1) * 2*pi/sum(QUANT1) )/2;
end

% PATCH
nAng = size(Angg,2);
Angg(2:nAng) = Angg(1:(nAng-1));
Angg(1) = 0;
Angg = Angg * 2*pi/Angg(end);

% another patch
Diss = [0, cumsum(vecnorm( diff(Bezz, 1, 2), 2, 1) )];

% report results
BezBase = Bezz;
AngBase = Angg;
DisBase = Diss;

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
% Display the control points of a path, and the orientation. 
%
% Special case for a hole in the sliding wheel. For convenience, the wheel
% is assumed to have radius 1 and will be properly scaled later.
%
% ---- INUPUT ------------------------------------------------------------
%  BPath         Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%
% ---- OUTPUT ------------------------------------------------------------
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function PlotHole( obj )

figure()
hold on
axis equal
grid on

ShowFig  = [obj.EvalAllPositions(), ...
  [cos(2*pi*(0:0.01:1)); sin(2*pi*(0:0.01:1))]];
fill(ShowFig(1,:),ShowFig(2,:), [.1,.1,.1], 'EdgeColor', 'none');

% plot control points
for i = 1:obj.nSegments
  CPts = obj.Segment{i}.CtrlPts;
  plot(CPts(1,[1,2]),CPts(2,[1,2]),'Color','blue')
  plot(CPts(1,[3,4]),CPts(2,[3,4]),'Color','blue')
  scatter(CPts(1,[2,3]), CPts(2,[2,3]),[],'blue','filled')
end

% points with numbers
for i = 1:2
  CPts = obj.Segment{i}.CtrlPts;
  scatter(CPts(1,1),CPts(2,1),200,'red','filled','o')
end

% center of wheel, for reference
scatter(0,0,200,'white','filled','o')

end
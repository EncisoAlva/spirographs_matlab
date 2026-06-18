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
function PlotHole( HPath )

figure()
hold on
axis equal
grid on

ShowFig  = [PathEval(HPath, 0.01), [cos(2*pi*(0:0.01:1)); sin(2*pi*(0:0.01:1))]];
fill(ShowFig(1,:),ShowFig(2,:), [.1,.1,.1], 'EdgeColor', 'none');

% plot control points
for i = 1:size(HPath,2)
  plot(HPath{i}(1,[1,2]),HPath{i}(2,[1,2]),'Color','blue')
  plot(HPath{i}(1,[3,4]),HPath{i}(2,[3,4]),'Color','blue')
  scatter(HPath{i}(1,[2,3]), HPath{i}(2,[2,3]),[],'blue','filled')
end

% points with numbers
for i = 1:2
  scatter(HPath{i}(1,1),HPath{i}(2,1),200,'red','filled','o')
  text(HPath{i}(1,1),HPath{i}(2,1),num2str(i-1),'Color','white','HorizontalAlignment','center')
end
for i = 3:size(HPath,2)
  scatter(HPath{i}(1,1),HPath{i}(2,1),200,'blue','filled','o')
  text(HPath{i}(1,1),HPath{i}(2,1),num2str(i-1),'Color','white','HorizontalAlignment','center')
end

% center of wheel, for reference
scatter(0,0,200,'white','filled','o')

end
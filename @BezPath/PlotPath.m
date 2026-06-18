% Display the control points of a path, and the orientation.
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
function PlotPath( BPath )

figure()
hold on
axis equal
grid on

ShowFig  = PathEval(BPath, 0.01);
fill(ShowFig(1,:),ShowFig(2,:), [.1,.1,.1], 'EdgeColor', 'none');

% plot control points
for i = 1:size(BPath,2)
  plot(BPath{i}(1,[1,2]),BPath{i}(2,[1,2]),'Color','blue')
  plot(BPath{i}(1,[3,4]),BPath{i}(2,[3,4]),'Color','blue')
  scatter(BPath{i}(1,[2,3]), BPath{i}(2,[2,3]),[],'blue','filled')
end

% normal vector
NormVec = EvalBezierNormal(BPath{1},0,0.5);
quiver(BPath{1}(1,1),BPath{1}(2,1),NormVec(1),NormVec(2),'Color','red')

% points with numbers
scatter(BPath{1}(1,1),BPath{1}(2,1),200,'red','filled','o')
text(BPath{1}(1,1),BPath{1}(2,1),'0','Color','white','HorizontalAlignment','center')
for i = 2:size(BPath,2)
  scatter(BPath{i}(1,1),BPath{i}(2,1),200,'blue','filled','o')
  text(BPath{i}(1,1),BPath{i}(2,1),num2str(i-1),'Color','white','HorizontalAlignment','center')
end

end
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
function PlotPath( obj )

figure()
hold on
axis equal
grid on

ShowFig  = obj.EvalPosition();
fill(ShowFig(1,:),ShowFig(2,:), [.1,.1,.1], 'EdgeColor', 'none');

% plot control points
for i = 1:obj.nSegments
  CurrSegment = obj.Segment{i};
  CPts = CurrSegment.CtrlPts;
  plot(CPts(1,[1,2]),CPts(2,[1,2]),'Color','blue')
  plot(CPts(1,[3,4]),CPts(2,[3,4]),'Color','blue')
  scatter(CPts(1,[2,3]), CPts(2,[2,3]),[],'blue','filled')
end

% normal vector
CurrSegment = obj.Segment{1};
CPts = CurrSegment.CtrlPts;
NormVec = CurrSegment.EvalNormal( 0,0.5 );
quiver(CPts(1,1),CPts(2,1),NormVec(1),NormVec(2),'Color','red')

% points with numbers
for i = [2:obj.nSegments, 1]
  CurrSegment = obj.Segment{i};
  CPts = CurrSegment.CtrlPts;
  scatter(CPts(1,1),CPts(2,1),200,'blue','filled','o')
  text(   CPts(1,1),CPts(2,1),num2str(i-1),'Color','white','HorizontalAlignment','center')
end
scatter(CPts(1,1),CPts(2,1),200,'red','filled','o')
text(   CPts(1,1),CPts(2,1),'0','Color','white','HorizontalAlignment','center')

end
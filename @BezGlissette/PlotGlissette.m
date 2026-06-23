function PlotGlissette( obj )

% process colors, due to updates
obj.ProcessColors()

% new figure
figure();
hold on
axis equal
axis on
grid on
set(gca,'Color','k')

% fill the bezier curve, for reference
fill(obj.DecorativeBez(1,:),obj.DecorativeBez(2,:), ...
  .15*[1,1,1], 'EdgeColor', 'none')

% actual plot
if obj.Multicolor
  drawPts = [ obj.MarkerPos, NaN(2,1)];
  drawCol = [ obj.ColorFunc'; NaN(1,3)];
  fill(drawPts(1,:),drawPts(2,:),...
    [0,0,0],'FaceVertexCData',drawCol,'EdgeColor','interp','LineWidth',2)
else
  plot(obj.MarkerPos(1,:),obj.MarkerPos(2,:),...
    'Color',obj.ColorVector{1}, 'LineWidth',2)
end

% reference points for multicolor
if obj.Multicolor
  switch obj.ColorRefCurve
    case 'CumDist'
      ColorVal = cumsum([0, vecnorm( diff( obj.MarkerPos, 1,2), 2,1 )]);
    case 'Bezier'
      ColorVal = cumsum([0, vecnorm( diff( obj.BezierPos, 1,2), 2,1 )]);
  end
  ColorVal = ColorVal/ColorVal(end);
  for rref = 0:(0.5/obj.ColorCycles):1
    [~,iidx] = max(ColorVal(ColorVal<=rref));
    scatter(obj.MarkerPos(1,iidx),obj.MarkerPos(2,iidx),'red','filled','o')
  end
end

end
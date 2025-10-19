%%
% fig 1: simple cycloid

% parameters
Tol = 0.001;

% line to roll on
LineL = { LineToBezier([0,0]',[2*pi,0]') };

% specific to this example
WheelRadius  = 1;
MarkerRadius = 1;

% Bezier curve
DecorativeBezier = [ [0-1,0]',[2*pi+1,0]' ];

% spirograph curve
[LocTime, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  AllBeziers_web( LineL, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve = zeros(1,size(BezierPos,2));
TimeFromCurve(2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
tick_x = 0:(2*pi/8):(2*pi);
tick_y = zeros(size(tick_x));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)];
aang = 0:(2*pi/8):(2*pi);

[~, i] = max(TimeFromCurve(TimeFromCurve<0.4));

% graph
figure(Theme="light")
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))];
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',1, 'Color','blue')
end
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')
%
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
scatter(tick_x, tick_y,'k', 'filled')

text([0,2*pi],[0,0]-.15,{'0','$2\pi R_\mathcal{W}$'},'Color','black','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+1.15,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')



addpath([pwd,'/webpage_animations'])

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
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( LineL, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
tick_x = 0:(2*pi/8):(2*pi);
tick_y = zeros(size(tick_x));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)];
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = [0:(1/120):1, 1-(0:(1/120):1)]

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))];
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text([0,2*pi],[0,0]-.15,{'0','$2\pi R_\mathcal{W}$'},'Color','black','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+1.15,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(2*pi+0.75,.15,'$\mathcal{L}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(pi,-.3,'Cycloid','Color','black','FontSize',20,'HorizontalAlignment','center', 'FontName','Helvetica')
scatter([-1,2*pi+1],[-.5,2.5],'white', 'filled')

text(0,max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')


xlim([-1,2*pi+1])
ylim([-1.7,3.7])

exportgraphics(f1,"fig01.gif",Append=true)

end

%%
% fig 7a: curtate cycloid

% parameters
Tol = 0.001;

% line to roll on
LineL = { LineToBezier([0,0]',[2*pi,0]') };

% specific to this example
WheelRadius  = 1;
MarkerRadius = 0.5;

% Bezier curve
DecorativeBezier = [ [0-1,0]',[2*pi+1,0]' ];

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( LineL, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
tick_x = 0:(2*pi/8):(2*pi);
tick_y = zeros(size(tick_x));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)];
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = [0:(1/120):1, 1-(0:(1/120):1)]

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))];
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text([0,2*pi],[0,0]-.15,{'0','$2\pi R_\mathcal{W}$'},'Color','black','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+1.15,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(2*pi+0.75,.15,'$\mathcal{L}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(pi,-.3,'Curtate cycloid','Color','black','FontSize',20,'HorizontalAlignment','center', 'FontName','Helvetica')
scatter([-1,2*pi+1],[-.5,2.5],'white', 'filled')

text(0,max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')

xlim([-1,2*pi+1])
ylim([-1.7,3.7])

exportgraphics(f1,"fig07a.gif",Append=true)

end

%%
% fig 7b: prolate cycloid

% parameters
Tol = 0.001;

% line to roll on
LineL = { LineToBezier([0,0]',[2*pi,0]') };

% specific to this example
WheelRadius  = 1;
MarkerRadius = 1.5;

% Bezier curve
DecorativeBezier = [ [0-1,0]',[2*pi+1,0]' ];

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( LineL, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
tick_x = 0:(2*pi/8):(2*pi);
tick_y = zeros(size(tick_x));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)];
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = [0:(1/120):1, 1-(0:(1/120):1)]

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))];
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text([0,2*pi],[0,0]-.15,{'0','$2\pi R_\mathcal{W}$'},'Color','black','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+1.15,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(2*pi+0.75,.15,'$\mathcal{L}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')

text(pi,-.3,'Prolate cycloid','Color','black','FontSize',20,'HorizontalAlignment','center', 'FontName','Helvetica')
scatter([-1,2*pi+1],[-.5,2.5],'white', 'filled')

xlim([-1,2*pi+1])
ylim([-1.7,3.7])

exportgraphics(f1,"fig07b.gif",Append=true)

end

%%
% figure 2a

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/2;
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/16):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/240):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 2$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3],...
  [-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig02a.gif",Append=true)

end

%%
% figure 2b

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/3;
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*3)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/240):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],...
  [-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig02b.gif",Append=true)

end

%%
% figure 2c

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/5;
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*5)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/240):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 5$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],...
  [-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig02c.gif",Append=true)

end

%%
% figure 3a

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(3/2);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(4*3)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(2*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 3/2$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],...
  [-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig03a.gif",Append=true)

end

%%
% figure 3b

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(5/2);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(4*5)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(2*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 5/2$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],...
  [-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig03b.gif",Append=true)

end

%%
% figure 3c

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(7/3);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 3);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*7)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(3*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1-.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-max(WheelRadius+MarkerRadius)-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 7/3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],...
  [-1-(WheelRadius+MarkerRadius)-.3, 1+(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig03c.gif",Append=true)

end

%%
% figure 5a

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
Circ = FlipPath(Circ);
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(3);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

BezierPos = flip(BezierPos,2);
WhCtrPos = flip(WhCtrPos,2);
MarkerPos = flip(MarkerPos,2);
MarkerAngle = flip(MarkerAngle,2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*3)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(1*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)-.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1+.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-.15-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-.3, 1+.3],...
  [-1-.3, 1+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig05a.gif",Append=true)

end

%%
% figure 5b

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
Circ = FlipPath(Circ);
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(5);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

BezierPos = flip(BezierPos,2);
WhCtrPos = flip(WhCtrPos,2);
MarkerPos = flip(MarkerPos,2);
MarkerAngle = flip(MarkerAngle,2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*5)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(1*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)-.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1+.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-.15-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 5$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-.3, 1+.3],...
  [-1-.3, 1+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig05b.gif",Append=true)

end

%%
% figure 5c

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
Circ = FlipPath(Circ);
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(5/3);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 3);

BezierPos = flip(BezierPos,2);
WhCtrPos = flip(WhCtrPos,2);
MarkerPos = flip(MarkerPos,2);
MarkerAngle = flip(MarkerAngle,2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*5)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(3*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)-.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1+.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-.15-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 5/3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

%xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
%ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-.3, 1+.3],...
  [-1-.3, 1+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig05c.gif",Append=true)

end

%%
% fig 5d

% parameters
Tol = 0.001;

% line to roll on
Circ = { [...
  [-1, 0]',...
  [-1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 1, 0]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 1, 0]'...
  ],[...
  [ 1, 0]',...
  [ 1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]'...
  [ 0,-1]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ] };
Circ = FlipPath(Circ);
%PlotBezierCtrlPts(Circ)

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(6/5);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 5);

BezierPos = flip(BezierPos,2);
WhCtrPos = flip(WhCtrPos,2);
MarkerPos = flip(MarkerPos,2);
MarkerAngle = flip(MarkerAngle,2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% line ticks
aang = 0:(2*pi/(8*3)):(2*pi);
tick_x = cos(aang);
tick_y = sin(aang);

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);


f1 = figure(Theme="light");
for t = 0:(1/(5*240)):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
%text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)-.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,1+.15,0,'$\mathcal{C}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,-1-.15-.15,0,'$R_\mathcal{C}/R_\mathcal{W} = 6/5$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])
ylim([-1-max(WheelRadius+MarkerRadius)-.3, 1+max(WheelRadius+MarkerRadius)+.3])

scatter([-1-.3, 1+.3],...
  [-1-.3, 1+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig05d.gif",Append=true)

end

%%
% figure 4a

% parameters
Tol = 0.001;

% line to roll on
Circ = struct2cell(load('ExampleCurves.mat','Bone'));
Circ = Circ{1};
%PlotBezierCtrlPts(Circ)

nCurves = size(Circ,2);

% shrink the curve to het a redius of approx 2pi
Perim = PathPerimeter(Circ,Tol);
for i = 1:nCurves
  Circ{i} = Circ{i}*(2*pi/Perim);
end

% fun with perimeters
Perim = zeros(1,size(Circ,2));
for i = 1:size(Circ,2)
  Perim(i) = PathPerimeter({Circ{i}},Tol);
end
PerimCum = [0,cumsum(Perim)];

Tvals = 0:(Tol):1;
CumArc = zeros(nCurves,size(Tvals,2));
for i = 1:nCurves
  CurrBez = EvalBezier(Circ{i}, Tvals);
  dists = vecnorm( diff(CurrBez,1,2), 2, 1);
  CumArc(i,2:end) = cumsum(dists);
end

% line ticks
aang = 0:(2*pi/(8*3)):(2*pi);
tick_x = zeros(1,size(aang,2));
tick_y = zeros(1,size(aang,2));
for i = 1:size(aang,2)
  [m,j]  = max(PerimCum(PerimCum<=aang(i)));
  j = mod(j-1,nCurves)+1;
  res = aang(i) - m;
  [~,j2] = max( CumArc(j, CumArc(j,:) <= res ) );
  t_approx = ( Tvals(j2)*(CumArc(j,j2+1)-res) + Tvals(j2+1)*(res-CumArc(j,j2)) )/( CumArc(j,j2+1)-CumArc(j,j2) );
  Ptmp = EvalBezier( Circ{j}, t_approx );
  tick_x(i) = Ptmp(1);
  tick_y(i) = Ptmp(2);
end

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/3;
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 1);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);

LowBorder = min(DecorativeBezier(2,:));
TopBorder = max(DecorativeBezier(2,:));
RRBorder = max(DecorativeBezier(1,:));
LLBorder = min(DecorativeBezier(1,:));

f1 = figure(Theme="light");
for t = 0:(1/480):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,TopBorder-.15,0,'$\mathcal{B}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,LowBorder-(WheelRadius+MarkerRadius)-.15,0,'$P_\mathcal{B}/P_\mathcal{W} = 3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3])
ylim([LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3])

scatter([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3],...
  [LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig04a.gif",Append=true)

end

%%
% figure 4b

% parameters
Tol = 0.001;

% line to roll on
Circ = struct2cell(load('ExampleCurves.mat','Bone'));
Circ = Circ{1};
%PlotBezierCtrlPts(Circ)

nCurves = size(Circ,2);

% shrink the curve to het a redius of approx 2pi
Perim = PathPerimeter(Circ,Tol);
for i = 1:nCurves
  Circ{i} = Circ{i}*(2*pi/Perim);
end

% fun with perimeters
Perim = zeros(1,size(Circ,2));
for i = 1:size(Circ,2)
  Perim(i) = PathPerimeter({Circ{i}},Tol);
end
PerimCum = [0,cumsum(Perim)];

Tvals = 0:(Tol):1;
CumArc = zeros(nCurves,size(Tvals,2));
for i = 1:nCurves
  CurrBez = EvalBezier(Circ{i}, Tvals);
  dists = vecnorm( diff(CurrBez,1,2), 2, 1);
  CumArc(i,2:end) = cumsum(dists);
end

% line ticks
aang = 0:(2*pi/(8*5/2)):(2*pi);
tick_x = zeros(1,size(aang,2));
tick_y = zeros(1,size(aang,2));
for i = 1:size(aang,2)
  [m,j]  = max(PerimCum(PerimCum<=aang(i)));
  j = mod(j-1,nCurves)+1;
  res = aang(i) - m;
  [~,j2] = max( CumArc(j, CumArc(j,:) <= res ) );
  t_approx = ( Tvals(j2)*(CumArc(j,j2+1)-res) + Tvals(j2+1)*(res-CumArc(j,j2)) )/( CumArc(j,j2+1)-CumArc(j,j2) );
  Ptmp = EvalBezier( Circ{j}, t_approx );
  tick_x(i) = Ptmp(1);
  tick_y(i) = Ptmp(2);
end

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(5/2);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 2);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);

LowBorder = min(DecorativeBezier(2,:));
TopBorder = max(DecorativeBezier(2,:));
RRBorder = max(DecorativeBezier(1,:));
LLBorder = min(DecorativeBezier(1,:));

f1 = figure(Theme="light");
for t = 0:(2/480):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,TopBorder-.15,0,'$\mathcal{B}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,LowBorder-(WheelRadius+MarkerRadius)-.15,0,'$P_\mathcal{B}/P_\mathcal{W} = 5/2$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3])
ylim([LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3])

scatter([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3],...
  [LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig04b.gif",Append=true)

end

%%
% figure 4c

% parameters
Tol = 0.001;

% line to roll on
Circ = struct2cell(load('ExampleCurves.mat','Bone'));
Circ = Circ{1};
%PlotBezierCtrlPts(Circ)

nCurves = size(Circ,2);

% shrink the curve to het a redius of approx 2pi
Perim = PathPerimeter(Circ,Tol);
for i = 1:nCurves
  Circ{i} = Circ{i}*(2*pi/Perim);
end

% fun with perimeters
Perim = zeros(1,size(Circ,2));
for i = 1:size(Circ,2)
  Perim(i) = PathPerimeter({Circ{i}},Tol);
end
PerimCum = [0,cumsum(Perim)];

Tvals = 0:(Tol):1;
CumArc = zeros(nCurves,size(Tvals,2));
for i = 1:nCurves
  CurrBez = EvalBezier(Circ{i}, Tvals);
  dists = vecnorm( diff(CurrBez,1,2), 2, 1);
  CumArc(i,2:end) = cumsum(dists);
end

% line ticks
aang = 0:(2*pi/(8*7)):(2*pi);
tick_x = zeros(1,size(aang,2));
tick_y = zeros(1,size(aang,2));
for i = 1:size(aang,2)
  [m,j]  = max(PerimCum(PerimCum<=aang(i)));
  j = mod(j-1,nCurves)+1;
  res = aang(i) - m;
  [~,j2] = max( CumArc(j, CumArc(j,:) <= res ) );
  if j2==size(Tvals,2)
    t_approx = 1;
  else
    t_approx = ( Tvals(j2)*(CumArc(j,j2+1)-res) + Tvals(j2+1)*(res-CumArc(j,j2)) )/( CumArc(j,j2+1)-CumArc(j,j2) );
  end
  Ptmp = EvalBezier( Circ{j}, t_approx );
  tick_x(i) = Ptmp(1);
  tick_y(i) = Ptmp(2);
end

% specific to this example
WheelRadius  = (PathPerimeter(Circ,Tol)/(2*pi))/(7/3);
MarkerRadius = WheelRadius;

% Bezier curve
DecorativeBezier = PathEval(Circ, Tol);

% spirograph curve
[~, BezierPos, WhCtrPos, MarkerPos, MarkerAngle] = ...
  GenerateGlissette_web( Circ, WheelRadius, MarkerRadius, pi, ...
    Tol, Tol, 3);

% make time based on the traveres arc length
TimeFromCurve_pre = zeros(2,size(BezierPos,2));
TimeFromCurve_pre(1,2:end) = cumsum( vecnorm( diff(BezierPos,1,2), 2, 1 ) );
TimeFromCurve_pre(2,2:end) = cumsum( vecnorm( diff(WhCtrPos, 1,2), 2, 1 ) );
TimeFromCurve = mean(TimeFromCurve_pre,1);
TimeFromCurve = TimeFromCurve*(1/TimeFromCurve(end));

% circle, the picture
aang = 0:(2*pi/50):(2*pi);
circle_pts = [cos(aang);sin(aang)]*WheelRadius;
aang = 0:(2*pi/8):(2*pi);

LowBorder = min(DecorativeBezier(2,:));
TopBorder = max(DecorativeBezier(2,:));
RRBorder = max(DecorativeBezier(1,:));
LLBorder = min(DecorativeBezier(1,:));

f1 = figure(Theme="light");
for t = 0:(1.5/480):1

[~, i] = max(TimeFromCurve(TimeFromCurve<=t));

% graph
clf(f1)
plot(DecorativeBezier(1,:), DecorativeBezier(2,:), 'LineWidth',2, 'Color','black')
axis equal
hold on
axis off
plot(MarkerPos(1,:), MarkerPos(2,:), 'LineWidth',2, 'Color','red')
text(0, 1+max(WheelRadius+MarkerRadius,2*WheelRadius)+.15,'JCEA2025','Color',[.95,.95,.95],'FontSize',16,'HorizontalAlignment','left', 'Interpreter','latex')
%
axle_pts = [cos(aang + MarkerAngle(i)); sin(aang + MarkerAngle(i))]*WheelRadius;
for j = 1:size(axle_pts,2)
  plot(WhCtrPos(1,i)+[0,axle_pts(1,j)],WhCtrPos(2,i)+[0,axle_pts(2,j)], 'LineWidth',.1, 'Color','blue')
end
plot([WhCtrPos(1,i),MarkerPos(1,i)],[WhCtrPos(2,i),MarkerPos(2,i)], 'LineWidth',1.5, 'Color','blue')
scatter(WhCtrPos(1,i), WhCtrPos(2,i),'blue', 'filled')
plot(circle_pts(1,:)+WhCtrPos(1,i), circle_pts(2,:)+WhCtrPos(2,i), 'LineWidth',1.5, 'Color','blue')
%
scatter(tick_x, tick_y,'k', 'filled')
%
scatter(BezierPos(1,i), BezierPos(2,i),'blue', 'filled')
scatter(MarkerPos(1,i), MarkerPos(2,i), 75,'red', 'filled')

text(MarkerPos(1,i),MarkerPos(2,i)+.15,'$\mathbf{M}$','Color','red','FontSize',16,'HorizontalAlignment','center', 'Interpreter','latex')
text(WhCtrPos(1,i)+0.15+WheelRadius,WhCtrPos(2,i),'$\mathcal{W}$','Color','blue','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')
text(0,TopBorder-.15,0,'$\mathcal{B}$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

text(0,LowBorder-(WheelRadius+MarkerRadius)-.15,0,'$P_\mathcal{B}/P_\mathcal{W} = 7/3$','Color','black','FontSize',18,'HorizontalAlignment','center', 'Interpreter','latex')

xlim([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3])
ylim([LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3])

scatter([ LLBorder-max(WheelRadius+MarkerRadius)-.3,  RRBorder+max(WheelRadius+MarkerRadius)+.3],...
  [LowBorder-max(WheelRadius+MarkerRadius)-.3, TopBorder+max(WheelRadius+MarkerRadius)+.3],[],[.95,.95,.95])

exportgraphics(f1,"fig04c.gif",Append=true)

end
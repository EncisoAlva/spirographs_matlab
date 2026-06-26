%%
% generate points according to the following algorithm:
% 1. create points uniformly on a square that coers the curve
% 2. consider only points inside the curve
% 3. compute the distance curve-point, if the point is far from the curve,
%    then it is more likely to be removed
% 4. comb the points based on its distance to the curve

% parameters
nStars = 5000;
starRadius = 0.001;

% create box containing the curve
Xmin =  Inf;
Xmax = -Inf;
Ymin =  Inf;
Ymax = -Inf;
for i = 1:nPts
  Xmin = min( Xmin, min(AllMarkerPos{i}(1,:)) );
  Xmax = max( Xmax, max(AllMarkerPos{i}(1,:)) );
  Ymin = min( Ymin, min(AllMarkerPos{i}(2,:)) );
  Ymax = max( Ymax, max(AllMarkerPos{i}(2,:)) );
end

% main loop
rng(2026)
WB = waitbar(0,strcat('Generating points...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2026)');
%
coordStar = zeros(2, nStars);
currStar = 1;
while currStar < nStars
  %
  % 
  if getappdata(WB,'canceling')
    disp('Ended by user.')
    close(v);
    delete(WB)
    break
  end
  %
  xp = Xmin + (Xmax - Xmin) * rand(1, 1);
  yp = Ymin + (Ymax - Ymin) * rand(1, 1);
  %
  % if point in curve
  if inpolygon(xp, yp, AllMarkerPos{1}(1,:), AllMarkerPos{1}(2,:))
    distStar = min( vecnorm( AllMarkerPos{1} - [xp; yp], 2, 1 ) );
    probKeep = 1 - distStar/(distStar+starRadius);
    %
    % trimming
    if rand(1, 1) <= probKeep
      coordStar(:,currStar) = [xp; yp];
      %
      % update progressbar
      waitbar(currStar/nStars,WB);
      currStar = currStar + 1;
    end
  end
end

% delete waitbar, if it is still around
if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% check
figure()
hold on
axis equal
grid on
set(gca,'color', 'k');
for i = 1:nPts
  plot(AllMarkerPos{i}(1,:),AllMarkerPos{i}(2,:),'Color',ColorVector{mod(i-1,k)+1}, 'LineWidth',2)
end
scatter(coordStar(1,:),coordStar(2,:),3,[48, 133, 194]/255,'filled','o')

% fancy options
x0 = min(DecorativeBez(1,:));
xF = max(DecorativeBez(1,:));
y0 = min(DecorativeBez(2,:));
yF = max(DecorativeBez(2,:));
for p = 1:nPts
  x0 = min( x0, min(AllMarkerPos{p}(1,:)) );
  xF = max( xF, max(AllMarkerPos{p}(1,:)) );
  y0 = min( y0, min(AllMarkerPos{p}(2,:)) );
  yF = max( yF, max(AllMarkerPos{p}(2,:)) );
end
%
x_ran = xF - x0;
y_ran = yF - y0;
%
if y_ran/x_ran < ExpectedRatio
  y_ran_new = x_ran*ExpectedRatio;
  y0 = y0 - (y_ran_new - y_ran)/2;
  yF = yF + (y_ran_new - y_ran)/2;
else 
  if y_ran/x_ran > ExpectedRatio
    x_ran_new = y_ran/ExpectedRatio;
    x0 = x0 - (x_ran_new - x_ran)/2;
    xF = xF + (x_ran_new - x_ran)/2;
  end
end
%
xlim([x0 xF])
ylim([y0 yF])

% the actual figure
figure()
hold on
axis equal
axis off
grid off
set(gca,'color', 'k');
scatter(coordStar(1,:),coordStar(2,:),2,[48, 133, 194]/255,'filled','o')


xlim([x0 xF])
ylim([y0 yF])

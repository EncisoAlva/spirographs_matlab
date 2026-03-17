% load curve
BPath_pack1 = struct2cell(load('ExampleCollections.mat','Circlegon'));
BPath_pack2 = BPath_pack1{1};
BPath0 = BPath_pack2{4};

clear BPath_pack1 BPath_pack2

BPath0 = ShiftPath( BPath0, 2, false );

PlotPath(BPath0)

%%
BPath1 = cell(1,16);

for i = 1:16
  th0 = 2*pi*(i-1)/16 + pi/2;
  th1 = 2*pi*(i  )/16 + pi/2;
  %
  P0 = 2*[cos(th0); sin(th0)];
  P3 = 2*[cos(th1); sin(th1)];
  %
  P1 = P0 + [cos(th0+pi/2); sin(th0+pi/2)]*(8/3)*tan(2*pi/(4*16));
  P2 = P3 + [cos(th1-pi/2); sin(th1-pi/2)]*(8/3)*tan(2*pi/(4*16));
  %
  BPath1{i} = [P0, P1, P2, P3];
end

PlotPath(BPath1)

clear i P0 P1 P2 P3 th0 th1

%%
DirVec0 = cell(1,16);
DirVec1 = cell(1,16);
DirMag0 = cell(1,16);
DirMag1 = cell(1,16);

for i = 1:16
  DirVec0{i} = [BPath0{i}(:,2) - BPath0{i}(:,1), BPath0{i}(:,3) - BPath0{i}(:,4)];
  DirVec1{i} = [BPath1{i}(:,2) - BPath1{i}(:,1), BPath1{i}(:,3) - BPath1{i}(:,4)];
  DirMag0{i} = [norm(DirVec0{i}(:,1)),norm(DirVec0{i}(:,2))];
  DirMag1{i} = [norm(DirVec1{i}(:,1)),norm(DirVec1{i}(:,2))];
end

%%
% parameters

% technical stuff
Tol = 0.005;
CloseTol = 0.01;
MaxSpins = 100;
WheelRadiusTol = 0.000001;

% designer stuff
MarkerAngle0 = 0;

WheelBezRatio = 12/5;

WheelMarkerRatio = 4/5;

% willing to loose 1% of total area due to each corner rounding
CornerRoundingRadius = sqrt(0.001*PathArea(BPath0, 0.005)/(pi));

%%
k = 0.2;

BPathK = cell(1,16);

for i = 1:16
  P0 = (1-k)*BPath0{i}(:,1) + k*BPath1{i}(:,1);
  P3 = (1-k)*BPath0{i}(:,4) + k*BPath1{i}(:,4);
  %
  DirMag01 = (1-k)*DirMag0{i}(1) + k*DirMag1{i}(1);
  DirMag34 = (1-k)*DirMag0{i}(2) + k*DirMag1{i}(2);
  %
  DirVec01 = (1-k)*DirVec0{i}(:,1) + k*DirVec1{i}(:,1);
  DirVec34 = (1-k)*DirVec0{i}(:,2) + k*DirVec1{i}(:,2);
  %
  DirVec01 = DirMag01 * DirVec01 / norm(DirVec01);
  DirVec34 = DirMag34 * DirVec34 / norm(DirVec34);
  %
  P1 = P0 + DirVec01;
  P2 = P3 + DirVec34;
  %
  BPathK{i} = [P0, P1, P2, P3];
end

PlotPath(BPathK)

if k~=1
  [BPath_rounded_flipped] = ...
    RemoveAllCorners( FlipPath(BPathK), CornerRoundingRadius, 0.005, false );
  BPathK = FlipPath(BPath_rounded_flipped);

  BPathK = ShiftPath( BPathK, 2, true );
end

PlotPath(BPathK)

%%
MaxTime = 4;

VidName = 'Morhping01';

ColorVector = { [255, 59, 209]/255 };

%%

fps     = 30;
MaxTime = ceil(MaxTime*fps)/fps;
nTimes  = MaxTime*fps;

CurveOpts = {};
CurveOpts.CloseEnds = false;
CurveOpts.Tol = Tol;
CurveOpts.CloseTol = CloseTol;

aang = 2*pi*(0:1/1:1) + pi;
aang(end) = [];
MarkerAngle0Array = aang;
nPts = size(MarkerAngle0Array,2);
k = size(ColorVector,2);

%%

% original figure
close all
set(0, 'DefaultFigureColor', 'k');

f1 = figure('Visible','off','Name','Just the curve');
hold on
axis equal
axis off

ExpectedRatio = 16/9;
%
x0 = -2.1;
xF =  2.1;
y0 = -2.1;
yF =  2.1;
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
if ExpectedRatio == 16/9
  set(f1,'PaperPosition',[0 0 [1080 1920]*4],'PaperUnits','points');
end
%
f2 = figure('Visible','off','Name','With circle');

v = VideoWriter(strcat(VidName,".mp4"),'MPEG-4');
v.Quality = 100;
open(v)

% main loop
WB = waitbar(0,strcat('Generating video (',VidName,'.mp4)...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2025)');
for i = 0:nTimes
  % 
  % update progressbar
  waitbar(i/nTimes,WB);
  if getappdata(WB,'canceling')
    disp('Ended by user.')
    close(v);
    delete(WB)
    break
  end
  %
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)

  %
  k = i/nTimes;

  BPathK = cell(1,16);

  for j = 1:16
    P0 = (1-k)*BPath0{j}(:,1) + k*BPath1{j}(:,1);
    P3 = (1-k)*BPath0{j}(:,4) + k*BPath1{j}(:,4);
    %
    DirMag01 = (1-k)*DirMag0{j}(1) + k*DirMag1{j}(1);
    DirMag34 = (1-k)*DirMag0{j}(2) + k*DirMag1{j}(2);
    %
    DirVec01 = (1-k)*DirVec0{j}(:,1) + k*DirVec1{j}(:,1);
    DirVec34 = (1-k)*DirVec0{j}(:,2) + k*DirVec1{j}(:,2);
    %
    DirVec01 = DirMag01 * DirVec01 / norm(DirVec01);
    DirVec34 = DirMag34 * DirVec34 / norm(DirVec34);
    %
    P1 = P0 + DirVec01;
    P2 = P3 + DirVec34;
    %
    BPathK{j} = [P0, P1, P2, P3];
  end

  %PlotPath(BPathK)

  %if k~=1
  %  [BPath_rounded_flipped] = ...
  %    RemoveAllCorners( FlipPath(BPathK), CornerRoundingRadius, 0.005, false );
  %  BPathK = FlipPath(BPath_rounded_flipped);
%
%    BPathK = ShiftPath( BPathK, 2, true );
%  else
    BPathK = ShiftPath( BPathK, 2, false );
%  end

  %PlotPath(BPathK)

  WheelRadius = (PathPerimeter(BPathK,0.001)/(2*pi))/WheelBezRatio;
  MarkerRadius = WheelRadius*WheelMarkerRatio;

  % compute curves
  [ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
      SetupCurves_Npts( nPts, BPathK, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
        5, CurveOpts);

  plot(DecorativeBez(1,:),DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)

  % add all strokes of the marker up to the current time
  for p = 1:nPts
    plot(AllMarkerPos{p}(1,:),AllMarkerPos{p}(2,:),'Color',ColorVector{p},'LineWidth',2)
  end

  %
  writeVideo(v,getframe)
end

%%
% stop for some time after everything is finished, WITH the wheel on
for stopper = 0:(fps*(1))
  writeVideo(v,getframe)
end

% main loop
WB = waitbar(0,strcat('Generating video (',VidName,'.mp4)...'), ...
  'Name','Spirograph over Bezier curves by Enciso-Alva (2025)');
for i = 0:nTimes
  % 
  % update progressbar
  waitbar(i/nTimes,WB);
  if getappdata(WB,'canceling')
    disp('Ended by user.')
    close(v);
    delete(WB)
    break
  end
  %
  % keep the base elements
  clf(f2)
  copyobj(f1.Children,f2)
  set(0,"CurrentFigure",f2)

  %
  k = 1-i/nTimes;

  BPathK = cell(1,16);

  for j = 1:16
    P0 = (1-k)*BPath0{j}(:,1) + k*BPath1{j}(:,1);
    P3 = (1-k)*BPath0{j}(:,4) + k*BPath1{j}(:,4);
    %
    DirMag01 = (1-k)*DirMag0{j}(1) + k*DirMag1{j}(1);
    DirMag34 = (1-k)*DirMag0{j}(2) + k*DirMag1{j}(2);
    %
    DirVec01 = (1-k)*DirVec0{j}(:,1) + k*DirVec1{j}(:,1);
    DirVec34 = (1-k)*DirVec0{j}(:,2) + k*DirVec1{j}(:,2);
    %
    DirVec01 = DirMag01 * DirVec01 / norm(DirVec01);
    DirVec34 = DirMag34 * DirVec34 / norm(DirVec34);
    %
    P1 = P0 + DirVec01;
    P2 = P3 + DirVec34;
    %
    BPathK{j} = [P0, P1, P2, P3];
  end

  %PlotPath(BPathK)

  %if k~=1
  %  [BPath_rounded_flipped] = ...
  %    RemoveAllCorners( FlipPath(BPathK), CornerRoundingRadius, 0.005, false );
  %  BPathK = FlipPath(BPath_rounded_flipped);
%
%    BPathK = ShiftPath( BPathK, 2, true );
%  else
    BPathK = ShiftPath( BPathK, 2, false );
%  end

  %PlotPath(BPathK)

  WheelRadius = (PathPerimeter(BPathK,0.001)/(2*pi))/WheelBezRatio;
  MarkerRadius = WheelRadius*WheelMarkerRatio;

  % compute curves
  [ DecorativeBez, ~, ~, ~, AllMarkerPos, ~ ] = ...
      SetupCurves_Npts( nPts, BPathK, WheelRadius, MarkerRadius, MarkerAngle0Array, ...
        5, CurveOpts);

  plot(DecorativeBez(1,:),DecorativeBez(2,:),'Color',[.4 .4 .4],'LineWidth',2)

  % add all strokes of the marker up to the current time
  for p = 1:nPts
    plot(AllMarkerPos{p}(1,:),AllMarkerPos{p}(2,:),'Color',ColorVector{p},'LineWidth',2)
  end

  %
  writeVideo(v,getframe)
end

%%
% stop for some time after everything is finished, WITH the wheel on
for stopper = 0:(fps*(1))
  writeVideo(v,getframe)
end

% delete waitbar, if it is still around
if exist('WB','var')
  if getappdata(WB,'canceling')
    return
  end
  delete(WB)
end

% finalize the video object and close the figure
close(v);
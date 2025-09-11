MaxDisDelta = 0.001;
WheelRadius = 0.5;

Curve1 = [-1,1; -1+(4/3)*tan(pi/8),1; 0,(4/3)*tan(pi/8); 0,0]';
Curve2 = [0,0; 0,(4/3)*tan(pi/8); 1-(4/3)*tan(pi/8),1; 1,1]';

% initial and finishing times
t0_1 = 0;
t0_2 = 0;
tF_1 = 1;
tF_2 = 1;

iter = 1;

while true
LocalTime1 = t0_1:((tF_1-t0_1)/(ceil((tF_1-t0_1)/MaxDisDelta)+iter)):tF_1;
LocalTime2 = t0_2:((tF_2-t0_2)/(ceil((tF_2-t0_2)/MaxDisDelta)+iter)):tF_2;

BezierPos1  = EvalBezier( Curve1, LocalTime1 );
BezierPos2  = EvalBezier( Curve2, LocalTime2 );
BezierNorm1 = EvalBezierNormal( Curve1, LocalTime1, WheelRadius );
BezierNorm2 = EvalBezierNormal( Curve2, LocalTime2, WheelRadius );

WhCtrPos1 = BezierPos1 + BezierNorm1;
WhCtrPos2 = BezierPos2 + BezierNorm2;

figure()
hold on
plot(BezierPos1(1,:),BezierPos1(2,:))
plot(BezierPos2(1,:),BezierPos2(2,:))
plot(WhCtrPos1(1,:),WhCtrPos1(2,:))
plot(WhCtrPos2(1,:),WhCtrPos2(2,:))
axis equal

curve1_dis = zeros(size(LocalTime1));
curve1_dis_idx = zeros(size(LocalTime1));
for i = 1:size(LocalTime1,2)
  [dist, idxx] = min( vecnorm( WhCtrPos1(:,i) - WhCtrPos2, 2, 1 ) );
  curve1_dis(i) = dist;
  curve1_dis_idx(i) = idxx;
end

[dist, idx1] = min(curve1_dis);
idx2 = curve1_dis_idx(idx1);

plot([BezierPos1(1,idx1),WhCtrPos1(1,idx1)],[BezierPos1(2,idx1),WhCtrPos1(2,idx1)])
plot([BezierPos2(1,idx2),WhCtrPos2(1,idx2)],[BezierPos2(2,idx2),WhCtrPos2(2,idx2)])

if dist < MaxDisDelta
  break
end

t0_1 = max(0,LocalTime1(idx1)-0.5/2^iter);
t0_2 = max(0,LocalTime2(idx2)-0.5/2^iter);
tF_1 = min(1,LocalTime1(idx1)+0.5/2^iter);
tF_2 = min(1,LocalTime1(idx1)+0.5/2^iter);

iter = iter+1;
end

CrossTime1 = LocalTime1(idx1);
CrossTime2 = LocalTime1(idx2);

t1 = CrossTime1;
t2 = CrossTime2;

% Casteljau division at the point of intersection
Curve1_updated = zeros(2,4);
Curve2_updated = zeros(2,4);

Curve1_updated(:,1) = Curve1(:,1);
Curve1_updated(:,2) = (1-t1)*Curve1(:,1) + t1*Curve1(:,2);
Curve1_updated(:,3) = (1-t1)^2*Curve1(:,1) + 2*(1-t1)*t1*Curve1(:,2) + t1^2*Curve1(:,3);
Curve1_updated(:,4) = (1-t1)^3*Curve1(:,1) + 3*(1-t1)^2*t1*Curve1(:,2) + 3*(1-t1)*t1^2*Curve1(:,3) + t1^3*Curve1(:,4);

Curve2_updated(:,4) = Curve2(:,4);
Curve2_updated(:,3) = (1-t2)*Curve2(:,3) + t2*Curve2(:,4);
Curve2_updated(:,2) = (1-t2)^2*Curve2(:,2) + 2*(1-t2)*t2*Curve2(:,3) + t2^2*Curve2(:,4);
Curve2_updated(:,1) = (1-t2)^3*Curve2(:,1) + 3*(1-t2)^2*t2*Curve2(:,2) + 3*(1-t2)*t2^2*Curve2(:,3) + t2^3*Curve2(:,4);

% getting normal vectors to use the 4/3 tan(th/4) approx of a circle
N1 = EvalBezierNormal(Curve1,t1,1);
N2 = EvalBezierNormal(Curve2,t2,1);
GapAngle = atan2(N2(2),N2(1)) - atan2(N1(2),N1(1));

Curve_gap = zeros(2,4);
Curve_gap(:,1) = EvalBezier(Curve1,t1);
Curve_gap(:,4) = EvalBezier(Curve2,t2);
Curve_gap(:,2) = Curve_gap(:,1) + [0,1;-1,0]*N1* (4/3)*tan(GapAngle/4)*WheelRadius;
Curve_gap(:,3) = Curve_gap(:,4) + [0,-1;1,0]*N2* (4/3)*tan(GapAngle/4)*WheelRadius;

% plot
LocalTime = 0:(1/(ceil(1/MaxDisDelta)+iter)):1;

BezierPos1_new  = EvalBezier( Curve1_updated, LocalTime );
BezierPos2_new  = EvalBezier( Curve2_updated, LocalTime );
BezierGap  = EvalBezier( Curve_gap, LocalTime );

figure()
hold on
plot(BezierPos1_new(1,:),BezierPos1_new(2,:))
plot(BezierPos2_new(1,:),BezierPos2_new(2,:))
plot(BezierGap(1,:),BezierGap(2,:))
plot(BezierPos1(1,:),BezierPos1(2,:))
plot(BezierPos2(1,:),BezierPos2(2,:))
plot([BezierPos1(1,idx1),WhCtrPos1(1,idx1)],[BezierPos1(2,idx1),WhCtrPos1(2,idx1)])
plot([BezierPos2(1,idx2),WhCtrPos2(1,idx2)],[BezierPos2(2,idx2),WhCtrPos2(2,idx2)])
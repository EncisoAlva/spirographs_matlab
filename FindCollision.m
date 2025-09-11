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

CrossTime1 = LocalTime1(idx1);
CrossTime2 = LocalTime1(idx2);
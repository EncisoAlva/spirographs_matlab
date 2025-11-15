% Small collection of curves that I have used to make animations. Feel free
% to use them for exploring.
%
% The difference with the file 'ExampleCurves' is that collections have
% figures with rotational symmetry parametrized by an integr; most of them
% are polygon-like shapes.

%%
% polygon
Polygon = cell(1,8);

% make one single cycloid

Side_base = {LineToBezier([0,1]',[0,-1]')};

nSegments  = size(Side_base,2);

for N = 2:8
% rotate and concatenate cycloids
Polygon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N+pi/2;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi); sin(th+pi) ]/tan(pi/N);
  %
  Side_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Side_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Polygon_N{i+curr_side*nSegments} = Side_rotated{i};
  end
end

Polygon{N} = Polygon_N;
end

clear N bar curr_side i nSegments Polygon_N ROT Side_base Side_rotated th

%%
% polygon with cycloids as sides
Cycloidgon = cell(1,7);

% make one single cycloid
R  = 1;
b  = 1;
dt = pi/2;
t = (0:dt:(2*pi));

xy     = [ R*t-pi - b*sin(t); R - b*cos(t) ];
xy_der = [ R      - b*cos(t);     b*sin(t) ];

Cycloid_base = cell(1,size(t,2)-1);
for i = 2:size(t,2)
  del_i = t(i) - t(i-1);
  CurrCurve = zeros(2,4);
  CurrCurve(:,1) = xy(:,i-1);
  CurrCurve(:,2) = xy(:,i-1) + (1/3)*xy_der(:,i-1)*del_i;
  CurrCurve(:,3) = xy(:,i)   - (1/3)*xy_der(:,i)  *del_i;
  CurrCurve(:,4) = xy(:,i);
  Cycloid_base{i-1} = CurrCurve;
end
Cycloid_base = FlipPath(Cycloid_base);

nSegments  = size(Cycloid_base,2);

for N = 2:8
% rotate and concatenate cycloids
Cycloidgon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*pi/tan(pi/N);
  %
  Cycloid_rotated = cell(1,nSegments);
  for i = 1:size(Cycloid_base, 2)
    Cycloid_rotated{i} = ROT * Cycloid_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Cycloidgon_N{i+curr_side*nSegments} = Cycloid_rotated{i};
  end
end

Cycloidgon{N} = Cycloidgon_N;
end

clear b bar curr_side CurrCurve Cycloid_base Cycloid_rotated Cycloidgon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with cycloids as sides, inside out
Cycloidgon_in = cell(1,7);

% make one single cycloid
R  = 1;
b  = 1;
dt = pi/2;
t = (0:dt:(2*pi));

xy     = [ R*t-pi - b*sin(t); R - b*cos(t) ];
xy_der = [ R      - b*cos(t);     b*sin(t) ];

Cycloid_base = cell(1,size(t,2)-1);
for i = 2:size(t,2)
  del_i = t(i) - t(i-1);
  CurrCurve = zeros(2,4);
  CurrCurve(:,1) = xy(:,i-1);
  CurrCurve(:,2) = xy(:,i-1) + (1/3)*xy_der(:,i-1)*del_i;
  CurrCurve(:,3) = xy(:,i)   - (1/3)*xy_der(:,i)  *del_i;
  CurrCurve(:,4) = xy(:,i);
  Cycloid_base{i-1} = CurrCurve;
end
Cycloid_base = FlipPath(Cycloid_base);

for i = 1:size(Cycloid_base,2)
  Cycloid_base{i} = [1,0;0,-1]*Cycloid_base{i};
end

nSegments  = size(Cycloid_base,2);

for N = 2:8
% rotate and concatenate cycloids
Cycloidgon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*pi/tan(pi/N);
  %
  Cycloid_rotated = cell(1,nSegments);
  for i = 1:size(Cycloid_base, 2)
    Cycloid_rotated{i} = ROT * Cycloid_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Cycloidgon_N{i+curr_side*nSegments} = Cycloid_rotated{i};
  end
end

Cycloidgon_in{N} = Cycloidgon_N;
end

clear b bar curr_side CurrCurve Cycloid_base Cycloid_rotated Cycloidgon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with semicircles as sides
Circlegon = cell(1,8);

% make one single semicircle
Semicircle_base = {[...
  [ 1, 0]',...
  [ 1, 0]' + [0,1]'*(4/3)*tan(pi/16),...
  [ 1, 1]'/sqrt(2) + [1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 1, 1]'/sqrt(2)...
  ],[...
  [ 1, 1]'/sqrt(2),...
  [ 1, 1]'/sqrt(2) + [-1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 0, 1]' + [ 1,0]'*(4/3)*tan(pi/16),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]' + [-1,0]'*(4/3)*tan(pi/16),...
  [-1, 1]'/sqrt(2) + [ 1, 1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 1]'/sqrt(2)...
  ],[...
  [-1, 1]'/sqrt(2),...
  [-1, 1]'/sqrt(2) + [-1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 0]' + [ 0, 1]'*(4/3)*tan(pi/16),...
  [-1, 0]'...
  ]};

nSegments  = size(Semicircle_base,2);

for N = 2:8
% rotate and concatenate cycloids
Circlegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_base, 2)
    Semicircle_rotated{i} = ROT * Semicircle_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Circlegon_N{i+curr_side*nSegments} = Semicircle_rotated{i};
  end
end

Circlegon{N} = Circlegon_N;
end

clear b bar curr_side CurrCurve Semicircle_base Semicircle_rotated Circlegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with semicircles as sides
Circlegon_in = cell(1,8);

% make one single semicircle
Semicircle_base = {[...
  [ 1, 0]',...
  [ 1, 0]' + [0,1]'*(4/3)*tan(pi/16),...
  [ 1, 1]'/sqrt(2) + [1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 1, 1]'/sqrt(2)...
  ],[...
  [ 1, 1]'/sqrt(2),...
  [ 1, 1]'/sqrt(2) + [-1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 0, 1]' + [ 1,0]'*(4/3)*tan(pi/16),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]' + [-1,0]'*(4/3)*tan(pi/16),...
  [-1, 1]'/sqrt(2) + [ 1, 1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 1]'/sqrt(2)...
  ],[...
  [-1, 1]'/sqrt(2),...
  [-1, 1]'/sqrt(2) + [-1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 0]' + [ 0, 1]'*(4/3)*tan(pi/16),...
  [-1, 0]'...
  ]};
for i = 1:size(Semicircle_base,2)
  Semicircle_base{i} = [1,0;0,-1]*Semicircle_base{i};
end

nSegments  = size(Semicircle_base,2);

for N = 2:8
% rotate and concatenate cycloids
Circlegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_base, 2)
    Semicircle_rotated{i} = ROT * Semicircle_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Circlegon_N{i+curr_side*nSegments} = Semicircle_rotated{i};
  end
end

Circlegon_in{N} = Circlegon_N;
end

clear b bar curr_side CurrCurve Semicircle_base Semicircle_rotated Circlegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with semicircles as sides
Crossed_Circlegon_2 = cell(1,15);

% make one single semicircle
Side_base = {[...
  [ 1, 0]',...
  [ 1, 0]' + [0,1]'*(4/3)*tan(pi/16),...
  [ 1, 1]'/sqrt(2) + [1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 1, 1]'/sqrt(2)...
  ],[...
  [ 1, 1]'/sqrt(2),...
  [ 1, 1]'/sqrt(2) + [-1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 0, 1]' + [ 1,0]'*(4/3)*tan(pi/16),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]' + [-1,0]'*(4/3)*tan(pi/16),...
  [-1, 1]'/sqrt(2) + [ 1, 1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 1]'/sqrt(2)...
  ],[...
  [-1, 1]'/sqrt(2),...
  [-1, 1]'/sqrt(2) + [-1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 0]' + [ 0, 1]'*(4/3)*tan(pi/16),...
  [-1, 0]'...
  ]};

nSegments  = size(Side_base,2);

for N = 3:2:15
% rotate and concatenate cycloids
Curve_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(2*pi/N);
  %
  Side_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Side_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Curve_N{i+curr_side*nSegments} = Side_rotated{i};
  end
end
%PlotPath(Curve_N)

Crossed_Circlegon_2{N} = Curve_N;
end

clear bar curr_side Curve_N i N ROT Side_base Side_rotated th nSegments

%%
% polygon with cycloids as sides
Epicycloid = cell(1,7);

for N = 1:8
% make one single cycloid
R  = 1;
r  = 1/N; 
Rr = (R+r)/r;
dt = 2*pi/(N*4);
t = (0:dt:(2*pi));

xy     = [ (R+r)*sin(t) -    r*sin(Rr*t) ; (R+r)*cos(t) -    r*cos(Rr*t) ];
xy_der = [ (R+r)*cos(t) - Rr*r*cos(Rr*t) ; -(R+r)*sin(t) + Rr*r*sin(Rr*t) ];

CurrCycloid = cell(1,size(t,2)-1);
for i = 2:size(t,2)
  del_i = t(i) - t(i-1);
  CurrCurve = zeros(2,4);
  CurrCurve(:,1) = xy(:,i-1);
  CurrCurve(:,2) = xy(:,i-1) + (1/3)*xy_der(:,i-1)*del_i;
  CurrCurve(:,3) = xy(:,i)   - (1/3)*xy_der(:,i)  *del_i;
  CurrCurve(:,4) = xy(:,i);
  CurrCycloid{i-1} = CurrCurve;
end
CurrCycloid = FlipPath(CurrCycloid);

Epicycloid{N} = CurrCycloid;
end

clear CurrCurve CurrCycloid del_i dt i N r R Rr t y xy xy_der

%%
% hypocycloid, direct from parametric function
Hypocycloid = cell(1,8);

for N = 2:8
% make one single cycloid
R  = 1;
r  = 1/N; 
Rr = (R-r)/r;
dt = 2*pi/(N*2);
t = (0:dt:(2*pi));

xy     = [ (R-r)*sin(t) -    r*sin(Rr*t) ;  (R-r)*cos(t) +    r*cos(Rr*t) ];
xy_der = [ (R-r)*cos(t) - Rr*r*cos(Rr*t) ; -(R-r)*sin(t) - Rr*r*sin(Rr*t) ];

CurrCycloid = cell(1,size(t,2)-1);
for i = 2:size(t,2)
  del_i = t(i) - t(i-1);
  CurrCurve = zeros(2,4);
  CurrCurve(:,1) = xy(:,i-1);
  CurrCurve(:,2) = xy(:,i-1) + (1/3)*xy_der(:,i-1)*del_i;
  CurrCurve(:,3) = xy(:,i)   - (1/3)*xy_der(:,i)  *del_i;
  CurrCurve(:,4) = xy(:,i);
  CurrCycloid{i-1} = CurrCurve;
end
CurrCycloid = FlipPath(CurrCycloid);
CurrCycloid = ShiftPath(CurrCycloid, 1,false);

Hypocycloid{N} = CurrCycloid;
end

clear CurrCurve CurrCycloid del_i dt i N r R Rr t y xy xy_der

%%
% circle crowned with lines outside
Target = cell(1,8);

for N = 2:8
% make one single semicircle
Target_N = cell(1,3*N+1);

for i = 0:(N-1)
  ang_i  = (i-1)*2*pi/N;
  ang_i1 =    i *2*pi/N;
  magn   = 1.5;
  Target_N{3*i+1} = [...
    [cos(ang_i ),sin(ang_i )]',...
    [cos(ang_i ),sin(ang_i )]'+[cos(ang_i +pi/2),sin(ang_i +pi/2)]'*(4/3)*tan(pi/(2*N)),...
    [cos(ang_i1),sin(ang_i1)]'+[cos(ang_i1-pi/2),sin(ang_i1-pi/2)]'*(4/3)*tan(pi/(2*N)),...
    [cos(ang_i1),sin(ang_i1)]'...
    ];
  Target_N{3*i+2} = LineToBezier([cos(ang_i1),sin(ang_i1)]',     [cos(ang_i1),sin(ang_i1)]'*magn);
  Target_N{3*i+3} = LineToBezier([cos(ang_i1),sin(ang_i1)]'*magn,[cos(ang_i1),sin(ang_i1)]');
end

% change starting point
[c1, c2] = HalfBezierSingle(Target_N{1});
Target_N{1}   = c2;
Target_N{end} = c1;
Target_N = FlipPath(Target_N);

Target{N} = Target_N;
end

clear ang_i ang_i1 c1 c2 i magn N Target_N

%%
% circle with lines inside
Target_in = cell(1,8);

for N = 2:8
% make one single semicircle
Target_N = cell(1,3*N+1);

for i = 0:(N-1)
  ang_i  = (i-1)*2*pi/N +pi/N +pi/2;
  ang_i1 =    i *2*pi/N +pi/N +pi/2;
  %magn   = (1-1/N)/(1+1/pi);
  magn = 0.5;
  Target_N{3*i+1} = [...
    [cos(ang_i ),sin(ang_i )]',...
    [cos(ang_i ),sin(ang_i )]'+[cos(ang_i +pi/2),sin(ang_i +pi/2)]'*(4/3)*tan(pi/(2*N)),...
    [cos(ang_i1),sin(ang_i1)]'+[cos(ang_i1-pi/2),sin(ang_i1-pi/2)]'*(4/3)*tan(pi/(2*N)),...
    [cos(ang_i1),sin(ang_i1)]'...
    ];
  Target_N{3*i+2} = LineToBezier([cos(ang_i1),sin(ang_i1)]',     [cos(ang_i1),sin(ang_i1)]'*magn);
  Target_N{3*i+3} = LineToBezier([cos(ang_i1),sin(ang_i1)]'*magn,[cos(ang_i1),sin(ang_i1)]');
end

% change starting point
[c1, c2] = HalfBezierSingle(Target_N{1});
Target_N{1}   = c2;
Target_N{end} = c1;
Target_N = FlipPath(Target_N);

Target_in{N} = Target_N;
end

clear ang_i ang_i1 c1 c2 i magn N Target_N

%%
% polygon with squares as sides
Squaregon = cell(1,8);

% make one single semicircle
Square_base = {...
  LineToBezier([1,0]',[ 1,2]'),...
  LineToBezier([1,2]',[-1,2]'),...
  LineToBezier([-1,2]',[-1,0]')...
  };

nSegments  = size(Square_base,2);

for N = 2:8
% rotate and concatenate cycloids
Squaregon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Square_rotated = cell(1,nSegments);
  for i = 1:size(Square_base, 2)
    Square_rotated{i} = ROT * Square_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Squaregon_N{i+curr_side*nSegments} = Square_rotated{i};
  end
end

Squaregon{N} = Squaregon_N;
end

clear b bar curr_side CurrCurve Square_base Square_rotated Squaregon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with squares as sides
Trianglegon = cell(1,8);

% make one single semicircle
Side_base = {...
  LineToBezier([1,0]',[0,sqrt(3)]'),...
  LineToBezier([0,sqrt(3)]',[-1,0]')...
  };

nSegments  = size(Side_base,2);

for N = 2:8
% rotate and concatenate cycloids
Trianglegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Triangle_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Triangle_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Trianglegon_N{i+curr_side*nSegments} = Triangle_rotated{i};
  end
end

Trianglegon{N} = Trianglegon_N;
end

clear b bar curr_side CurrCurve Side_base Triangle_rotated Trianglegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with spikes as sides
Spikegon = cell(1,8);

% make one single semicircle
Side_base = {[...
  [ 1, 0]',...
  [ 1, 0]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0, 1]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]'+[ 0,-1]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ]};

nSegments  = size(Side_base,2);

for N = 2:8
% rotate and concatenate cycloids
Spikegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Side_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Side_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Spikegon_N{i+curr_side*nSegments} = Side_rotated{i};
  end
end

Spikegon{N} = Spikegon_N;
end

clear b bar curr_side CurrCurve Side_base Side_rotated Spikegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with spikes as sides
Spikegon_in = cell(1,8);

% make one single semicircle
Side_base = {[...
  [ 1, 0]',...
  [ 1, 0]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]',...
  [ 0,-1]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ]};

nSegments  = size(Side_base,2);

for N = 2:8
% rotate and concatenate cycloids
Spikegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]/tan(pi/N);
  %
  Side_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Side_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Spikegon_N{i+curr_side*nSegments} = Side_rotated{i};
  end
end

Spikegon_in{N} = Spikegon_N;
end

clear b bar curr_side CurrCurve Side_base Side_rotated Spikegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with spikes as sides
Spikegon_in2 = cell(1,8);

% make one single semicircle
Side_base = {...
  LineToBezier([2,0]',[1,0]'),...
  [...
  [ 1, 0]',...
  [ 1, 0]'+[-1, 0]'*(4/3)*tan(pi/8),...
  [ 0,-1]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]',...
  [ 0,-1]'+[ 0, 1]'*(4/3)*tan(pi/8),...
  [-1, 0]'+[ 1, 0]'*(4/3)*tan(pi/8),...
  [-1, 0]'...
  ],...
  LineToBezier([-1,0]',[-2,0]')...
  };

nSegments  = size(Side_base,2);

for N = 2:8
% rotate and concatenate cycloids
Spikegon_N = cell(1,nSegments*N);

for curr_side = 0:(N-1)
  th  = curr_side*2*pi/N;
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*2/tan(pi/N);
  %
  Side_rotated = cell(1,nSegments);
  for i = 1:size(Side_base, 2)
    Side_rotated{i} = ROT * Side_base{i} + bar;
  end
  % 
  for i = 1:nSegments
    Spikegon_N{i+curr_side*nSegments} = Side_rotated{i};
  end
end

Spikegon_in2{N} = Spikegon_N;
end

clear b bar curr_side CurrCurve Side_base Side_rotated Spikegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with semicircles as sides
TruncatedCirclegon = cell(1,8);

% make one single semicircle
Semicircle_big = {[...
  [ 2, 0]',...
  [ 2, 0]' + [0,1]'*(8/3)*tan(pi/16),...
  [ 2, 2]'/sqrt(2) + [1,-1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [ 2, 2]'/sqrt(2)...
  ],[...
  [ 2, 2]'/sqrt(2),...
  [ 2, 2]'/sqrt(2) + [-1,1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [ 0, 2]' + [ 1,0]'*(8/3)*tan(pi/16),...
  [ 0, 2]'...
  ],[...
  [ 0, 2]',...
  [ 0, 2]' + [-1,0]'*(8/3)*tan(pi/16),...
  [-2, 2]'/sqrt(2) + [ 1, 1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [-2, 2]'/sqrt(2)...
  ],[...
  [-2, 2]'/sqrt(2),...
  [-2, 2]'/sqrt(2) + [-1,-1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [-2, 0]' + [ 0, 1]'*(8/3)*tan(pi/16),...
  [-2, 0]'...
  ]};

Semicircle_small = {[...
  [ 1, 0]',...
  [ 1, 0]' + [0,1]'*(4/3)*tan(pi/16),...
  [ 1, 1]'/sqrt(2) + [1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 1, 1]'/sqrt(2)...
  ],[...
  [ 1, 1]'/sqrt(2),...
  [ 1, 1]'/sqrt(2) + [-1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 0, 1]' + [ 1,0]'*(4/3)*tan(pi/16),...
  [ 0, 1]'...
  ],[...
  [ 0, 1]',...
  [ 0, 1]' + [-1,0]'*(4/3)*tan(pi/16),...
  [-1, 1]'/sqrt(2) + [ 1, 1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 1]'/sqrt(2)...
  ],[...
  [-1, 1]'/sqrt(2),...
  [-1, 1]'/sqrt(2) + [-1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1, 0]' + [ 0, 1]'*(4/3)*tan(pi/16),...
  [-1, 0]'...
  ]};

nSegments  = size(Semicircle_big,2);

for N = 2:8
% rotate and concatenate cycloids
TruncatedCirclegon_N = cell(1,nSegments*2*N);

bar_short_length = ( 2+ 1/sin( pi*(N-2)/(2*N) ) ) / sin(pi/N) - cos(pi*(N-2)/(2*N))/sin(pi*(N-2)/(2*N));
bar_long_length  = ( 2+ 1/sin( pi*(N-2)/(2*N) ) ) * cos(pi/N) / sin(pi/N);

for curr_side = 0:(N-1)
  %
  % big circle
  th  = (2*curr_side)*2*pi/(2*N);
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*bar_long_length;
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_big, 2)
    Semicircle_rotated{i} = ROT * Semicircle_big{i} + bar;
  end
  % 
  for i = 1:nSegments
    TruncatedCirclegon_N{i+(2*curr_side)*nSegments} = Semicircle_rotated{i};
  end
  %
  % small circle
  th  = (2*curr_side+1)*2*pi/(2*N);
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*bar_short_length;
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_small, 2)
    Semicircle_rotated{i} = ROT * Semicircle_small{i} + bar;
  end
  % 
  for i = 1:nSegments
    TruncatedCirclegon_N{i+(2*curr_side+1)*nSegments} = Semicircle_rotated{i};
  end
end

TruncatedCirclegon{N} = ShiftPath(TruncatedCirclegon_N, 2, false);
end

clear b bar bar_short_length bar_long_length curr_side CurrCurve Semicircle_big Semicircle_small Semicircle_rotated TruncatedCirclegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with semicircles as sides
TruncatedCirclegon_alt = cell(1,8);

% make one single semicircle
Semicircle_big = {[...
  [ 2, 0]',...
  [ 2, 0]' + [0,1]'*(8/3)*tan(pi/16),...
  [ 2, 2]'/sqrt(2) + [1,-1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [ 2, 2]'/sqrt(2)...
  ],[...
  [ 2, 2]'/sqrt(2),...
  [ 2, 2]'/sqrt(2) + [-1,1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [ 0, 2]' + [ 1,0]'*(8/3)*tan(pi/16),...
  [ 0, 2]'...
  ],[...
  [ 0, 2]',...
  [ 0, 2]' + [-1,0]'*(8/3)*tan(pi/16),...
  [-2, 2]'/sqrt(2) + [ 1, 1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [-2, 2]'/sqrt(2)...
  ],[...
  [-2, 2]'/sqrt(2),...
  [-2, 2]'/sqrt(2) + [-1,-1]'*(8/3)*tan(pi/16)/sqrt(2),...
  [-2, 0]' + [ 0, 1]'*(8/3)*tan(pi/16),...
  [-2, 0]'...
  ]};

Semicircle_small = {[...
  [ 1, 0]',...
  [ 1, 0]' + [0,-1]'*(4/3)*tan(pi/16),...
  [ 1,-1]'/sqrt(2) + [1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 1,-1]'/sqrt(2)...
  ],[...
  [ 1,-1]'/sqrt(2),...
  [ 1,-1]'/sqrt(2) + [-1,-1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [ 0,-1]' + [ 1,0]'*(4/3)*tan(pi/16),...
  [ 0,-1]'...
  ],[...
  [ 0,-1]',...
  [ 0,-1]' + [-1,0]'*(4/3)*tan(pi/16),...
  [-1,-1]'/sqrt(2) + [ 1, -1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1,-1]'/sqrt(2)...
  ],[...
  [-1,-1]'/sqrt(2),...
  [-1,-1]'/sqrt(2) + [-1,1]'*(4/3)*tan(pi/16)/sqrt(2),...
  [-1,0]' + [ 0, -1]'*(4/3)*tan(pi/16),...
  [-1, 0]'...
  ]};

nSegments  = size(Semicircle_big,2);

for N = 2:8
% rotate and concatenate cycloids
TruncatedCirclegon_N = cell(1,nSegments*2*N);

bar_short_length = ( 2+ 1/sin( pi*(N-2)/(2*N) ) ) / sin(pi/N) - cos(pi*(N-2)/(2*N))/sin(pi*(N-2)/(2*N));
bar_long_length  = ( 2+ 1/sin( pi*(N-2)/(2*N) ) ) * cos(pi/N) / sin(pi/N);

for curr_side = 0:(N-1)
  %
  % big circle
  th  = (2*curr_side)*2*pi/(2*N);
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*bar_long_length;
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_big, 2)
    Semicircle_rotated{i} = ROT * Semicircle_big{i} + bar;
  end
  % 
  for i = 1:nSegments
    TruncatedCirclegon_N{i+(2*curr_side)*nSegments} = Semicircle_rotated{i};
  end
  %
  % small circle
  th  = (2*curr_side+1)*2*pi/(2*N);
  ROT = [cos(th), -sin(th); sin(th), cos(th)];
  bar = [ cos(th+pi/2); sin(th+pi/2) ]*bar_short_length;
  %
  Semicircle_rotated = cell(1,nSegments);
  for i = 1:size(Semicircle_small, 2)
    Semicircle_rotated{i} = ROT * Semicircle_small{i} + bar;
  end
  % 
  for i = 1:nSegments
    TruncatedCirclegon_N{i+(2*curr_side+1)*nSegments} = Semicircle_rotated{i};
  end
end

TruncatedCirclegon_alt{N} = ShiftPath(TruncatedCirclegon_N, 2, false);
end

clear b bar bar_short_length bar_long_length curr_side CurrCurve Semicircle_big Semicircle_small Semicircle_rotated TruncatedCirclegon_N del_i dt i nSegments R ROT t th Tol xy xy_der N

%%
% polygon with an unit circle pn each corner
Ballgon = cell(1,8);

for N = 2:8
  Section1 = [...
    [cos(-pi/N),sin(-pi/N)]'/tan(pi/N),...
    [cos(-pi/N),sin(-pi/N)]'*(1/tan(pi/N) + (4/3)*tan(pi*(N+2)/(8*N))),...
    [1,0]'*(1/sin(pi/N) + 1) + [0,-1]'*(4/3)*tan(pi*(N+2)/(8*N)),...
    [1,0]'*(1/sin(pi/N) + 1)...
    ];
  Section2 = [...
    [1,0]'*(1/sin(pi/N) + 1)...
    [1,0]'*(1/sin(pi/N) + 1) + [0,1]'*(4/3)*tan(pi*(N+2)/(8*N)),...
    [cos(pi/N),sin(pi/N)]'*(1/tan(pi/N) + (4/3)*tan(pi*(N+2)/(8*N))),...
    [cos(pi/N),sin(pi/N)]'/tan(pi/N)...
    ];

  Figure_N = cell(1,N*2);
  for i = 1:N
    th  = (i-1)*2*pi/N + pi/2;
    ROT = [cos(th), -sin(th); sin(th), cos(th)];
    Figure_N{2*i-1} = ROT * Section1;
    Figure_N{2*i  } = ROT * Section2;
  end

  Ballgon{N} = Figure_N;
end

clear Figure_N i N ROT Section1 Section2 th

%%
% self-intersecting polygon with an unit circle pn each corner
Angel = cell(1,16);

for N = 2:8
  Section1 = [...
    [cos(-pi/N),sin(-pi/N)]'/tan(pi/N),...
    [cos(-pi/N),sin(-pi/N)]'*(1/tan(pi/N) + (4/3)*tan(pi*(N+2)/(8*N))),...
    [1,0]'*(1/sin(pi/N) + 1) + [0,-1]'*(4/3)*tan(pi*(N+2)/(8*N)),...
    [1,0]'*(1/sin(pi/N) + 1)...
    ];
  Section2 = [...
    [1,0]'*(1/sin(pi/N) + 1)...
    [1,0]'*(1/sin(pi/N) + 1) + [0,1]'*(4/3)*tan(pi*(N+2)/(8*N)),...
    [cos(pi/N),sin(pi/N)]'*(1/tan(pi/N) + (4/3)*tan(pi*(N+2)/(8*N))),...
    [cos(pi/N),sin(pi/N)]'/tan(pi/N)...
    ];
  Section_in = [...
    [cos(-pi/N),sin(-pi/N)]'/tan(pi/N),...
    [cos(-pi/N),sin(-pi/N)]'*(1/tan(pi/N) - (4/3)*tan(pi*(N-2)/(4*N))),...
    [cos( pi/N),sin( pi/N)]'*(1/tan(pi/N) - (4/3)*tan(pi*(N-2)/(4*N))),...
    [cos( pi/N),sin( pi/N)]'/tan(pi/N)...
    ];

  if ceil(N/2)*2 - N > 0
    topp = 2*N;
    nCurves = 3*N;
  else
    topp = N;
    nCurves = 3*N/2;
  end

  Figure_N = cell(1,nCurves);
  counter = 1;
  for i = 1:topp
    th  = (i-1)*2*pi/N + pi/2;
    ROT = [cos(th), -sin(th); sin(th), cos(th)];
    if ceil(i/2)*2 - i > 0
      Figure_N{counter  } = ROT * Section1;
      Figure_N{counter+1} = ROT * Section2;
      counter = counter+2;
    else
      Figure_N{counter} = ROT * Section_in;
      counter = counter+1;
    end
  end

  Angel{N} = Figure_N;
end

clear counter Figure_N i N nCurves topp ROT Section1 Section2 Section_in th

%%
save('ExampleCollections.mat')

clear all
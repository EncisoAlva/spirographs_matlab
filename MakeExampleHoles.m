% Small collection of curves that I have used to make animations. Feel free
% to use them for exploring.
%
% This curves are intended to be used as holes inside a rolling wheel.
% Atlough several curves can be used, this file contain curves specifically
% designed for that purpose.

%%
% ice cream cone, semicircle facing up
IcecreamUp = cell(1,16);

for N = 2:16
  Figure_N = cell(1,4);
  Figure_N{1} = LineToBezier([0,0]',[ sin(pi/N),cos(pi/N)]');
  Figure_N{2} = [...
    [ sin(pi/N),cos(pi/N)]', ...
    [ sin(pi/N),cos(pi/N)]' + [-cos(pi/N),sin(pi/N)]'*tan(pi/(4*N)),...
    [0,1]' + [ 1,0]'*tan(pi/(4*N)),...
    [0,1]',...
    ];
  Figure_N{3} = [...
    [0,1]',...
    [0,1]' + [-1,0]'*tan(pi/(4*N)),...
    [-sin(pi/N),cos(pi/N)]' + [ cos(pi/N),sin(pi/N)]'*tan(pi/(4*N)),...
    [-sin(pi/N),cos(pi/N)]',...
    ];
  Figure_N{4} = LineToBezier([-sin(pi/N),cos(pi/N)]',[0,0]');

  Figure_N = ShiftPath( Figure_N, 2, false );

  IcecreamUp{N} = Figure_N;
end

clear Figure_N N

%%
% ice cream cone, semicircle facing down
IcecreamDown = cell(1,16);

for N = 2:16
  Figure_N = cell(1,4);
  Figure_N{1} = LineToBezier([0,0]',[ sin(pi/N),cos(pi/N)]');
  Figure_N{2} = [...
    [ sin(pi/N),cos(pi/N)]', ...
    [ sin(pi/N),cos(pi/N)]' + [-cos(pi/N),-sin(pi/N)]'*tan(pi/(4*N)),...
    [0,2*cos(pi/N)-1]' + [ 1,0]'*tan(pi/(4*N)),...
    [0,2*cos(pi/N)-1]',...
    ];
  Figure_N{3} = [...
    [0,2*cos(pi/N)-1]',...
    [0,2*cos(pi/N)-1]' + [-1,0]'*tan(pi/(4*N)),...
    [-sin(pi/N),cos(pi/N)]' + [ cos(pi/N),-sin(pi/N)]'*tan(pi/(4*N)),...
    [-sin(pi/N),cos(pi/N)]',...
    ];
  Figure_N{4} = LineToBezier([-sin(pi/N),cos(pi/N)]',[0,0]');

  Figure_N = ShiftPath( Figure_N, 2, false );

  IcecreamDown{N} = Figure_N;
end

clear Figure_N N

%%
save('ExampleHoles.mat')

clear all
% Display the control points of a path, and the orientation.
%
% ---- INUPUT ------------------------------------------------------------
%  CtrlPtsArray  Array with control points for each one of the Bezier
%                curves that make the curve {?} <- [2,4]'s
%
% ---- OUTPUT ------------------------------------------------------------
%  CtrlPtsArrayFlipped  Array with control points for each one of the 
%                Bezier IN REVERSE ORDER {?} <- [2,4]'s
%
% The last control point of the last curve must be equal to the first
% control point of the first curve. This is not checked.
%
function PlotBezierCtrlPts( CtrlPtsArray )

figure()
hold on
axis equal
grid on

for i = 1:size(CtrlPtsArray,2)
  scatter(CtrlPtsArray{i}(1,:), CtrlPtsArray{i}(2,:))
end

scatter(CtrlPtsArray{1}(1,1),CtrlPtsArray{1}(2,1),'red','filled','o')
NormVec = EvalBezierNormal(CtrlPtsArray{1},0,1);
plot(CtrlPtsArray{1}(1,1)+[0,NormVec(1)],CtrlPtsArray{1}(2,1)+[0,NormVec(2)],'red')

end
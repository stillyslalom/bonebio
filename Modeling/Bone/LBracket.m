function [geom,Aeq,beq,Aineq,bineq,LB,UB] = LBracket(params)
% cutouts within a rectangle of size L & H
L = 1.0;%length
H = 0.4;%height
h = 0.1; % for loading 
r = 0.05;
tol = L/100;
if (nargin == 0)
    params = [0.1 0.1 0.1 0.1];
end
params = params(:); % column
a = params(1); % left corner cutout
b = params(2); % right corner cutout
c = params(3); % width of top cutout
d = params(4); % depth of top cutout
% constraints of the form Aeq*params = beq, and Aineq*params <= bineq
Aeq = []; beq = [];
Aineq = [];
bineq = [];
LB = tol*ones(4,1);

UB(1) = H*0.9;
UB(2) = H*0.9;
UB(3) = H*0.9;
UB(4) = 0.9*(L-H);

% now create the brep structure
geom.vertices = [a 0; L-b 0; L b;L H; L-h H; H+r H;H H+r; H L;H/2+c/2 L; H/2+c/2 L-d; H/2-c/2 L-d; H/2-c/2 L;0 L; 0 a; H+r H+r]';
geom.segments = zeros(4,14);
geom.segments(1,:) = 1; % all linear
geom.segments(1,6) = 2; % except this
geom.segments(2,:) = 1:14;
geom.segments(3,:) = [2:14,1];
geom.segments(4,6) = 15; % center
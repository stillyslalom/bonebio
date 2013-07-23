function [geom,Aeq,beq,Aineq,bineq,LB,UB] = cantileverBeam(params)
% cutouts within a rectangle of size L & H
L = 2.0;%length
H = 1;%height
h = 0.05; % loading 
tol = L/100;
if (nargin == 0)
    params = [0.1; 0.1];  
end
params = params(:); % column
a = params(1); % corner cutouts
b = params(2); % left edge cutout
% constraints of the form Aeq*params = beq, and Aineq*params <= bineq
Aeq = []; beq = [];
Aineq = [];
bineq = [];
LB(1) = tol;UB(1) = H/2-h/2-tol;
LB(2) = tol;UB(2) = H/2-tol;
if ((a >= H/2-h/2-tol) || (b >= H/2-tol))
    geom = [];
    return;
end

% now create the brep structure
geom.vertices = [b 0;0 -b; 0 -H/2;L-a -H/2;L -H/2+a; L -h/2; L h/2; L H/2-a; L-a H/2; 0 H/2; 0 b ]';
geom.segments = [1 1 2 0;1 2 3 0;1 3 4 0; 1 4 5 0; 1 5 6 0; 1 6 7 0; 1 7 8 0; 1 8 9 0; 1 9 10 0; 1 10 11 0; 1 11 1 0]';
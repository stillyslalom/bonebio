clc; clear all;  format compact; format long
problem = 3;
nElements = 500;
shape = 'Quadratic';
C_OBuo=0.01;% original concentration of uncommitted osteoblasts
C_OBpo=0.001;% original concentration of osteoblast precursors
C_OBao=0.0005;% original concentration of active osteoblasts 
C_OCpo=0.001;% original concentation of osteoclast precursors
C_OCao=0.0001;% original conentration of active osteoclasts

Phio = 0.05;% Cortical bone porosity
elem = 1:nElements;
a_P_OBp = 0.1; %"preosteoblastic proliferation fraction"
D_Pivonka_OBu = 7e-2; % maximum differentiation rate of osteoblast progenitor cells by Pivonka et al., (d^-1)
A_OBa = 2.1107e-1; %maximum apoptosis rate of active osteoblasts
D_OBp = 1.657e-1; %maximum differentiation rate of active osteoblasts
K_res = 2; %bone resorption rate 
Min_Pi_MECH_act = .5;
lambda = 4;


% K_PTH_REP_OB=2.226e-1;% repressor equlibruim constant (pM)
% K_PTH_ACT_OB=1.5e2;% activation equlibruim constant(pM)
% D_PTH=8.6e1;% degredation rate (1/d)
% B_PTH=2.5e2;% intrinsic production rate (pM/d)
% P_PTH_d=5e4;% PTH dosage (pM)

K_TGF_REP_OBp=1.7543e-4;% repressor equlibruim constant (pM)
K_TGF_ACT_OBu=5.6328e-4;% activation equlibruim constant(pM)
D_TGF_beta = 1;% degredation rate (1/d)
P_TGF_d=5e4;% TGF dosage (pM)

b=Bonebio5_7_13('cantileverBeam.brep',nElements,shape,'PlaneStress');
% b=b.setPTH(P_PTH_d,B_PTH,D_PTH,K_PTH_ACT_OB,K_PTH_REP_OB);
b=b.setCellConc(C_OBuo,C_OBpo,C_OBao,C_OCpo,C_OCao);
b=b.setRates(D_OBp,A_OBa,K_res,Min_Pi_MECH_act,lambda);
b=b.setTGFact(a_P_OBp,D_Pivonka_OBu,D_TGF_beta,K_TGF_ACT_OBu,K_TGF_REP_OBp,P_TGF_d);
b=b.setPorosity(Phio);
b=b.setPseudoDensity(b.myNumElems,1);
b=b.fixEdge(6);
b=b.fixEdge(7);
%b=b.fixEdge(1);
b=b.applyYForceOnEdge(3,1);
%b=b.solveFEProblem(); 
b=b.simulateBone();

% figure(4)
% b.plotMesh();
% b=b.solveFEProblem();
% figure(5)
% b.plotStress();
% figure(6);
 %b.plotTGF()

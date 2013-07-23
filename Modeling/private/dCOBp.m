function rate = dCOBp(DOBu,COBu,piActOBu,POBp,COBp,PIActOBp,DOBp,piRepOBp)
% [DOBu,COBu,piActOBu,POBp,COBp,PIActOBp,DOBp,piRepOBp] = argvec

% rate = av(1)*av(2)*av(3)+av(4)*av(5)*av(6)-av(7)*av(2)*av(8);
rate = DOBu*COBu*piActOBu+POBp*COBp*PIActOBp-DOBp*COBu*piRepOBp;

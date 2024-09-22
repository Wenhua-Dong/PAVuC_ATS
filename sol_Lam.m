function [Lam] = sol_Lam(Q,X,A,G0)

v = length(X);

Lam = cell(1,v);
for i = 1:v
    E = Q{i}*X{i}-A*G0{i};
    EF = sqrt(sum(E.*E, 1)) + eps;
    Lam{i} = diag(0.5./EF);
end
end
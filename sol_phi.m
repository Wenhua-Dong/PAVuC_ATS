function [phi] = sol_phi(Q,X,A,G0,alpha)

v = length(X);
tmp = 1/(1-alpha);
err = zeros(v,1);

for i=1:v    
    E = Q{i}*X{i} - A*G0{i};
    E21 = sum(sqrt(sum(E.*E, 1))) + eps;
    err(i) = E21^tmp;
end
phi = err./sum(err);
end

function [A] = sol_A(Q,X,G0,Lam,pa)

v = length(X);
Asum = 0;
for i = 1:v
    Asum = Asum + pa(i)*G0{i}*Lam{i}*X{i}'*Q{i}';
end
[Ua,~,Va] = svd(Asum,'econ');
At = Ua*Va';
A = At';
end
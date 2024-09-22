function [Q] = sol_Q(X,A,G0,Lam)

v = length(X);
Q = cell(1,v);
for i=1:v
    XLGtAt = X{i}*Lam{i}*G0{i}'*A';
    [Uq,~,Vq] = svd(XLGtAt,'econ');
    Qt = Uq*Vq';
    Q{i} = Qt';
end
end
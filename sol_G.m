function [G] = sol_G(Q,X,Lam,A,G0,Pi,na,pa,idt,mu)

v = length(X);
n = size(X{1},2);
m = size(A,2);
G = cell(1,v);

Gsum = 0;
for i = 1:v
    if i~=idt
        Gsum = Gsum+G0{i}(:,Pi{i});
    end
end

for i = 1:v
    if i~=idt
        % Calculate the transpose of Pi and convert it into an equivalent vector
        Pimtv = Pi{i};        
        Pim = zeros(n-na);         
        for j = 1:n-na
            Pim(Pi{i}(na+j)-na,j) = 1;
        end        
        [~,id] = max(Pim,[],2);
        Pimtv(na+1:n) = id+na;
        
        tmp = (pa(i)*A'*Q{i}*X{i}*Lam{i}+mu*G0{idt}(:,Pimtv))./repmat(pa(i)*diag(Lam{i})'+mu,m,1);
    else        
        tmp = (pa(i)*A'*Q{i}*X{i}*Lam{i}+mu*Gsum)./repmat(pa(i)*diag(Lam{i})'+mu*(v-1),m,1);
    end
    for h = 1:n
        ut = tmp(:,h);
        G{i}(:,h) = EProjSimplex_new(ut');
    end
end
end
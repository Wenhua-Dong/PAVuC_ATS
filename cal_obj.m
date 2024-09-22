function [term] = cal_obj(Q,X,A,G0,Pi,pa,idt,mu)

v = length(X);

term = 0;
for i = 1:v
    E = Q{i}*X{i}-A*G0{i};
    term = term + pa(i)*sum(sqrt(sum(E.*E, 1)));
    if i~=idt
        term = term + mu*norm(G0{i}(:,Pi{i})-G0{idt},'fro')^2;
    end
end
end
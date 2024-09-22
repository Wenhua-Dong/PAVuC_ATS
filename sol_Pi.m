function [Pi] = sol_Pi(G0,na,idt,ids)

v = length(G0);
n = size(G0{1},2);
Pi = repmat({(1:n)'},1,v);
S = cell(1,v);
for i = 1:v
    S{i} = G0{i}(:,na+1:end);
end

for i = 1:v
    if i~=idt
        % Calculate the 2-step probability, where
        % the initial state probability is removed.
        score = S{i}'*S{idt};
        for j = 1:n-na
            [~,idi] = max(score(ids{i}(j),:));
            Pi{i}(idi+na) = ids{i}(j)+na;
            score(:,idi) = -1;
        end
    end
end
end
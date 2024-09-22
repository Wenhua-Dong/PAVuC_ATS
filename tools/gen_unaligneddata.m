function [data,ind] = gen_unaligneddata(X,na)

v = length(X);
n = size(X{1},2);
data = cell(1,v);
ind = cell(1,v);

oind = (1:n);                    % The indices of the original data
rand('twister',5489);
aind = randperm(n,na);           % The indices of the aligned samples
oind(aind) = [];                 % The remaining indices
% Generate the shuffled data
for i = 1:v
    sind = randperm(n-na); 
    ind{i} = [aind, oind(sind)]; % Generate the indices of the shuffled data
    data{i} = X{i}(:,ind{i});
end
end
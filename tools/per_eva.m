function [mresult] = per_eva(G,Pi,k,gnd)

v = length(G);
Gfusion = zeros(size(G{1}));
for i = 1:v
    Gfusion = Gfusion+G{i}(:,Pi{i})/v;
end
[U,~,~]=svd(Gfusion','econ');
F = U(:,1:k);
F = F./repmat(sqrt(sum(F.^2,2)),1,size(F,2));

stream = RandStream.getGlobalStream;
reset(stream);
MAXiter = 200; % Maximum number of iterations for k-means
REPlic = 20;   % Number of replications for KMeans
result = zeros(5,3);
for j = 1 : 5
    label = kmeans(F, k, 'maxiter', MAXiter, 'replicates', REPlic, 'emptyaction', 'singleton');
    result(j, : ) = measurement(label,gnd); 
end
mresult = mean(result);
mstd = std(result);
end
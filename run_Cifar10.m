clear; close all; clc;
addpath(genpath('./'));
datadir='./datasets/';
load('Cifar10');

v = length(X);             % The number of views
n = length(Y);             % The number of samples
k = length(unique(Y));     % The number of clusters
rho = 0;                   % Alignment ratio
na = round(rho*n);         % The number of aligned samples

%% Parameter setting
m = 1*k;                   % The number of anchors
mu = 0.1;                  % Trade-off parameter
alpha = 1.1;               % Control parameter

%% Normalization
for i=1:v
    X{i} = zscore(X{i});
    X{i} = X{i}';
end

%% Generate the unaligned data
if na<n
    [data,ind] = gen_unaligneddata(X,na);
else
    data = X;
end
clear X

%% Run PAVuC-ATS
tic;
[G,Pi,idt] = PAVuC_ATS(data,k,m,na,mu,alpha); 
time = toc;

%% Performance evaluation
if na<n
    gnd = Y(ind{idt});
else
    gnd = Y;
end
result = per_eva(G,Pi,k,gnd);
fprintf('\n Results: ACC: %.4f, NMI: %.4f, F: %.4f', result(1), result(2), result(3));


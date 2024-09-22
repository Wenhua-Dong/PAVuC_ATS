function [G0,Pi,idt] = PAVuC_ATS(X,k,m,na,mu,alpha)

% ~~~~~~~~~~~~~~~~~~~~~~~~ Objective function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% sum_{i=1}^{v} (phi_{i})^{alpha}*||Q_{i}X_{i}-AG_{i}||_{2,1}
%             + mu*sum_{i\neq t}||G_{i}Pi_{i}-G_{t}||_{F}^{2},
%
% Q_{i}:  Projection matrix;
% X_{i}:  The i-th view. Each column is an observation;
% A:      Consistent anchors;
% G_{i}:  Anchor graph;
% Pi_{i}：Permutation;
% phi_{i} Weight factor;
% When verifying the experimental results in the manuscript, please remove
% the random number generator (i.e., line 9: rand('twister',5489)) from the
% gen_unaligneddata function.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

v = length(X);
n = size(X{1},2);
dl = k;            % Dimension of the latent space
maxIter = 60 ;     % The number of iterations
%% Initialization
A = rand(dl,m); 
G0 = repmat({rand(m,n)},1,v); 
Pi = repmat({(1:n)'},1,v);
Lam = repmat({eye(n)},1,v);
phi = ones(1,v)./v;
pa = phi.^alpha;
idt = 1;           % The index of the alignment template

ids = cell(1,v);   % The alignment order of the latent representations
for i = 1:v
    tmp = var(X{i}(:,na+1:n),0,1);
    tmp = tmp./sum(tmp);
    [~,ids{i}] = sort(tmp,"descend");
end

%% Run PAVuC-ATS
for iter = 1:maxIter
    %% Solve Q    
    [Q] = sol_Q(X,A,G0,Lam);

    %% Solve A
    [A] = sol_A(Q,X,G0,Lam,pa);

    %% Solve G
    [G] = sol_G(Q,X,Lam,A,G0,Pi,na,pa,idt,mu); 
    G0 = G;

    %% Solve Lambda
    [Lam] = sol_Lam(Q,X,A,G0);

    %% Solve Pi
    if na<n
        [Pi] = sol_Pi(G0,na,idt,ids);
    end

    %% Solve phi
    [phi] = sol_phi(Q,X,A,G0,alpha);

    %% Calculate the objective value
    [obj(iter)] = cal_obj(Q,X,A,G0,Pi,pa,idt,mu);

    if (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-7)
        break;
    end
    pa = phi.^alpha;
    [~,idv] = sort(phi,"descend");
    idt = idv(1);
end
end


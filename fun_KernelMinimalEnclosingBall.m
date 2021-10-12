function [outputArg1,outputArg2,outputArg3] = fun_KernelMinimalEnclosingBall(X,sigma,varargin)
%fun_KernelMinimalEnclosingBall Inputs a dataset NxD and outputs the
%support vectors as well as function handle that is true if evaluated
%inside the given set and false otherwise.

if nargin<3
    verbosity=1;
else
    verbosity=varargin{1};
end

%% Compute mutual Euclidean distance between the vectors (for Gaussian kernel)
X_X2=pdist2(X,X,'squaredeuclidean'); % fast computation of matrix of squared eucledian distances

%% Build the Gaussian kernel

% sigma=1;
K = exp(-X_X2/sigma^2); % Gaussian kernel
% K = (X*X'+1)^10;
% K=X*X';
%%  Prepare for the optimization 

H = 2*K; % Hessian
f = -diag(K); % Linear term

[N,D] = size(X);

% Equality constraints Aeq*x = beq
Aeq = zeros(N);
Aeq(1,:) = ones(1,N);
beq = zeros(N,1);
beq(1,1) = 1;

% Lower bound x>0
lb = zeros(N,1);

%% Solve the quadratic programming

tic
[alp,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,[]);

if verbosity>1
    sprintf('Quadprog time: %f',toc)
end

%% Obtain the support vectors
tol = 1e-4;
supports = find(abs(alp)>tol);

%% Function handles for the decision boundary

k_x = @(x) exp(-sum((X-x)'.^2)'/sigma^2);
% k_xx= @(x,y) exp(-gamma*sum((x-y)'.^2)); 
f_x = @(x) 1 - 2*k_x(x)'*alp - alp'*f + 2 * alp'*K*alp;
i_x = @(x) f_x(X(supports(1),:)) > f_x(x);

%% Outputs
outputArg1 = supports;
outputArg2 = f_x;
outputArg3 = i_x;
end


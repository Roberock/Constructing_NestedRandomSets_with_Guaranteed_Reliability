%% Generate data
DGM_type=3;
%%
% KernelType= 'rbf';
% Sigma=[1 0;0 1];
% N=30;
% % rng default %  for reproducibility
% XXX = mvnrnd([0,0], Sigma, N);
% X = XXX;
% if DGM_type==1 % parabola data
%     X(:,2)=X(:,1).^2+X(:,2);
% elseif DGM_type==2 %  Two moons data
%     X(N/2+1:end,2)=-((X(N/2+1:end,1)/2).^2+X(N/2+1:end,2))+10; %  Two moons data
% elseif DGM_type==3 % Signature ?
%     X(:,2)=X(:,1).^2+X(:,2);
%     X(N/2+1:end,2)=-(X(N/2+1:end,1)) ;
% elseif DGM_type==4 % simple
%     X = abs(X);
% end
%%
load('X')
%%

% H = zeros(N,N);
% for i=1:N
%     for j=1:N
%         H(i,j)=X(i,:)*X(j,:)';
%     end
% end
%X_X = zeros(N,N);
%for i = 1:N
%    for j = 1:N
%        X_X(i,j)=(X(i,:) - X(j,:))*(X(i,:) - X(j,:))';
%    end
%end
X_X=pdist2(X,X).^2; % fast computation of matrix of squared eucledian distances
%%

gamma=0.01;
K = exp(-gamma*X_X);
% K = (X*X'+1)^10;
% K=X*X';


H = 2*K;
f = -diag(K);

Aeq = zeros(N);
Aeq(1,:) = ones(1,N);
beq = zeros(N,1);
beq(1,1) = 1;

lb = zeros(N,1);

%%

x0 = zeros(N,1);
[alp,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,[]);

%%

supports = find(abs(alp)>1e-4);
%%
center = alp'*X;
R2 = alp'*diag((center-X)*(center-X)');

%%
p = nsidedpoly(1000, 'Center', center, 'Radius', sqrt(R2));
figure,
% plot(p, 'FaceColor', 'r')
hold on
scatter(X(:,1),X(:,2)), 
hold on, 
scatter(X(supports,1),X(supports,2),'r')
% for i=1:length(supports)
%  scatter(X(supports(i),1),X(supports(i),2),'r')
% end
axis equal
%% KERNEL
R2 =  1/2 * alp'*H*alp - alp'*f
%% Decision boundary

k_x = @(x) exp(-gamma*sum((X-x)'.^2)');
k_xx= @(x,y) exp(-gamma*sum((x-y)'.^2)); 

f_x = @(x) 1 - 2*k_x(x)'*alp - alp'*f + 2 * alp'*K*alp;


%%

MESH = linspace(-2.5,2.5);
[MX,MY] = meshgrid(MESH);

Z = zeros(100);
for i =1:100
    for j = 1:100
       Z(i,j) = f_x([MX(i,j),MY(i,j)]);
    end
end

%%
figure
contour(MX,MY,Z,[f_x(X(supports(3),:)),f_x(X(supports(3),:))])
hold on 
hold on
scatter(X(:,1),X(:,2)), 
hold on, 
scatter(X(supports,1),X(supports,2),'r') 
%for i=1:length(supports)
%scatter(X(supports(i),1),X(supports(i),2),'r')
%end
axis equal



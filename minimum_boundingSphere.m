function [R,C,LAMBDA]=minimum_boundingSphere(delta,Alpha)
% delta = the samples of the uncertain quantities [Nsample X Ndimensins]

N=size(delta,1);
Nd=size(delta,2);
%Alpha=0.9; % percent of retained samples (1-alpha) discarded
Confun= @(x) SphereConstraint(x,delta,N,Nd,Alpha); 
Objfun= @(x) x(1);
x0=[max(sum((delta-mean(delta)).^2,2)), mean(delta)];
lb=[0, -inf*ones(1,Nd)];
ub=+inf*ones(1,Nd+1);
[x_opt,~,~,~,LAMBDA]  = fmincon(Objfun,x0,[],[],[],[],lb,ub,Confun);
R=x_opt(1);
C=x_opt(2:end);
end
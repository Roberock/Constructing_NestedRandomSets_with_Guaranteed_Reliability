%% Try SVM-One Class
clc; clear variables; close all

%% Generate data
DGM_type=3;
KernelType= 'rbf';
Sigma=[8 0.4;0.4 1];
N=5000;
rng default %  for reproducibility
X = mvnrnd([0,0], Sigma, N);
if DGM_type==1 % parabola data
    X(:,2)=X(:,1).^2+X(:,2);
elseif DGM_type==2 %  Two moons data
    X(N/2+1:end,2)=-((X(N/2+1:end,1)/2).^2+X(N/2+1:end,2))+10; %  Two moons data
elseif DGM_type==3 % Signature ?
    X(:,2)=X(:,1).^2+X(:,2);
    X(N/2+1:end,2)=-(X(N/2+1:end,1)) ;
elseif DGM_type==4 % banknote autentification data
    addpath([pwd '\Datasets_all' ])
    load('Banknote_autentication.mat')
    X=[Banknote_autentication.X(:,3:4)];
    N=size(X,1);
elseif DGM_type==5 % banknote autentification data    
    addpath('D:\PAPERS SUBMITTED and TO SUB\EngAplicationsAI\Extension to Journal\Matlab codes\ARAMIS challenge simplified')
    Temp=load('Training_Data.mat');
    X=Temp.TrainingData.Signals';
    X =X(Temp.TrainingData.Labels==1,[7,10]);
    N=size(X,1);
end

Y = ones(1,size(X,1));
[PseudoRadius,PseudoCenter,~]=minimum_boundingSphere(X,1);
% PseudoCenter=mean(X);
% PseudoRadius=max(sqrt(sum((X-PseudoCenter).^2,2)).*1.0);
X_vtmp= randsphere(size(X,1),size(X,2),PseudoRadius.*1.1,'sphere');
X_virtual=X_vtmp+PseudoCenter;
Y_virtual = -1*ones(1,size(X_virtual,1));
 Y_fake2class=[Y ,Y_virtual];
 X_fake2class=[X ; X_virtual];

% Y_fake2class=[Y ];
% X_fake2class=[X ];

figure(1)
gplotmatrix(X_fake2class,[],Y_fake2class)
%% train model
 CostFP=1;
 CostFN=100;
%'Cost',[0 CostFP; CostFN 0]
N_tmp=N;
% prepare for plots
if  size(X,2)==2
    Nxmesh=500;
    Nymesh=500;
    x = linspace(min(X(:,1))-5,max(X(:,1))+5,Nxmesh);
    y = linspace(min(X(:,2))-5,max(X(:,2))+5,Nymesh);
    [XX,YY] = meshgrid(x,y);
    X_mesh=[XX(:),YY(:)];
    EpsMatUpper=zeros(Nxmesh,Nymesh);
end
%%
Niter=15;
for i=1:Niter
    SVMModel = fitcsvm(X_fake2class,Y_fake2class,'KernelFunction',KernelType,'KernelScale','auto',...
        'BoxConstraint',1,'Verbose',1);
    
    Nsv(i)=sum(SVMModel.IsSupportVector(1:N_tmp));
    
    Beta=10^-4;
    [EpsL(i),EpsU(i)]=epsLU_fast(sum(Nsv(1:i)),N,Beta);
    
    % X2=X;
    % Y2=Y;
    % X2(SVMModel.IsSupportVector,:)=[];
    % Y2(SVMModel.IsSupportVector)=[];
    % SVMModel2 = fitcsvm(X2,Y2,'KernelFunction','rbf' ,'Verbose',1);
    % Nsv2=sum(SVMModel2.IsSupportVector)+Nsv1;
    %% plot
    if  size(X,2)==2 % if 2-d we can plot
        [Y_pred,Distance]=SVMModel.predict(X_mesh);
        Distance=Distance(:,1);
        
        Y_pred_save{i}=Y_pred;
        EpsMatUpper=max(EpsMatUpper,reshape((Y_pred>0).*(EpsU(i)),Nxmesh,Nymesh));
        figure(2)
        subplot(1,2,1)
        contour(x,y,reshape(Y_pred,Nxmesh,Nymesh));
        hold on
        if i==1
            scatter(X(:,1),X(:,2),'.k')
        end
        scatter(X(SVMModel.IsSupportVector(1:N_tmp),1),X( SVMModel.IsSupportVector(1:N_tmp),2),'+r')
        %scatter(X2(SVMModel2.IsSupportVector==1,1),X2(SVMModel2.IsSupportVector==1,2),'xb')
    end
    
    X_fake2class(SVMModel.IsSupportVector(1:N_tmp),:)=[];
    Y_fake2class(SVMModel.IsSupportVector(1:N_tmp))=[];
    N_tmp=N-Nsv(i);
end

if  size(X,2)==2
    figure(2)
    subplot(1,2,2)
    contour3(x,y,EpsMatUpper,'ShowText','on')
    xlabel('x_1')
    ylabel('x_2')
    zlabel('1-Belief')
end


%% uncertainty propagation
% the second-order response surface model of steel plate structures 
% Zhao Y, Deng Z, Zhang X. A robust stochastic model updating method with resampling processing[J]. Mechanical Systems and Signal Processing, 2020, 136(C).
f1 = @(x) 1.31+0.2152*x(:,1)-0.01455*x(:,2)-0.00002823*x(:,1).^2-0.0004878*x(:,1).*x(:,2)+0.0006576*x(:,2).^2;
f2 = @(x) 44.08+0.5145*x(:,1)-0.1333*x(:,2)-0.00002943*x(:,1).^2-0.004055*x(:,1).*x(:,2)+0.005588*x(:,2).^2;
 Model_Sim= @(x) [f1(x),f2(x)];
% 
% f1 = @(x) sin(x(:,1)).*sin(x(:,2));
% f2 = @(x) exp(abs(100 - sqrt(x(:,1).^2+x(:,2).^2)/pi));
% Model_Sim = @(x) [f2(x), -0.0001 * (abs(f1(x).*f2(x))+1).^0.1];




% Model_Simuulator= @(x) [  exp( abs(x(:,1)-x(:,2))) ,f2(x)];
for i=1:Niter 
X_in=X_mesh(Y_pred_save{i}(:)==1,:);
Out{i,:}= Model_Sim(X_in); 
figure(3)
scatter3(Out{i,:}(:,1),Out{i,:}(:,2),ones(1,size(Out{i,:},1)).*EpsU(i),'o','filled','MarkerFaceAlpha',0.3)
hold on
end
xlabel('f1')
ylabel('f2')

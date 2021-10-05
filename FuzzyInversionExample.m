%% Try SVM-One Class
clc; clear variables; close all

%% Generate data
DGM_type=1;
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
make_plots = false;

i = 1;
alpha(i) = 0;
while alpha(i) < 0.95

    SVMModel = fitcsvm(X_fake2class,Y_fake2class,'KernelFunction',KernelType,'KernelScale','auto',...
        'BoxConstraint',1,'Verbose',1);
    
    Nsv(i)=sum(SVMModel.IsSupportVector(1:N_tmp));
    
    Beta=10^-4;
    [EpsL(i),EpsU(i)]=epsLU_fast(sum(Nsv(1:i)),N,Beta);
    
    SVM_save{i} = SVMModel;
    alpha(i+1) = EpsU(i);

    % X2=X;
    % Y2=Y;
    % X2(SVMModel.IsSupportVector,:)=[];
    % Y2(SVMModel.IsSupportVector)=[];
    % SVMModel2 = fitcsvm(X2,Y2,'KernelFunction','rbf' ,'Verbose',1);
    % Nsv2=sum(SVMModel2.IsSupportVector)+Nsv1;
    %% plot
    if  size(X,2)==2 && make_plots % if 2-d we can plot
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
    i = i +1
end

alpha = EpsU;

if  size(X,2)==2 && make_plots
    figure(2)
    subplot(1,2,2)
    contour3(x,y,EpsMatUpper,'ShowText','on')
    xlabel('x_1')
    ylabel('x_2')
    zlabel('1-Belief')
end


%% Inverse problem. Started data was Gassian. Is it inside the inverted fuzzy structure?
% No sampling example



XRanges = [-15, 15];
YRanges = [-22, 25];

Nxmesh=500;
Nymesh=500;
x = linspace(XRanges(1), XRanges(2), Nxmesh);
y = linspace(YRanges(1), YRanges(2), Nymesh);
[XX,YY] = meshgrid(x,y);
X_mesh=[XX(:),YY(:)];

Membership_values = zeros(length(X_mesh),1);

XX = X_mesh;

if DGM_type==1 % parabola data
    XX(:,2)=XX(:,1).^2+XX(:,2);
elseif DGM_type==2 %  Two moons data
    XX(N/2+1:end,2)=-((XX(N/2+1:end,1)/2).^2+XX(N/2+1:end,2))+10; %  Two moons data
elseif DGM_type==3 % Signature ?
    XX(:,2)=XX(:,1).^2+XX(:,2);
    XX(N/2+1:end,2)=-(XX(N/2+1:end,1));
end

for i = 1:length(SVM_save)
    [inside,~] = SVM_save{i}.predict(XX);
    Membership_values(inside == 1) = alpha(i);
end

figure()
contourf(x,y,reshape(Membership_values, Nxmesh, Nymesh),'ShowText','on')

figure()
mesh(x,y,reshape(Membership_values,Nxmesh, Nymesh))

% 
% %% Inverse problem. Started data was Gassian. Is it inside the inverted fuzzy structure?
% % Sampling example
% 
% Nsamples_fuzzy = 10^4;
% ranges = [10, -10];
% 
% % inputs samples
% 
% Membership_values = zeros(Nsamples_fuzzy,1);
% 
% input_samples = rand(Nsamples_fuzzy, 2) * (ranges(2) - ranges(1)) + ranges(1);
% 
% XX = input_samples;
% 
% if DGM_type==1 % parabola data
%     XX(:,2)=XX(:,1).^2+XX(:,2);
% elseif DGM_type==2 %  Two moons data
%     XX(N/2+1:end,2)=-((XX(N/2+1:end,1)/2).^2+XX(N/2+1:end,2))+10; %  Two moons data
% elseif DGM_type==3 % Signature ?
%     XX(:,2)=XX(:,1).^2+XX(:,2);
%     XX(N/2+1:end,2)=-(XX(N/2+1:end,1));
% end
% 
% for i = 1:length(SVM_save)
%     [inside,~] = SVM_save{i}.predict(XX);
%     Membership_values(inside == 1) = alpha(i);
% end
% 
% figure()
% scatter3(XX(:,1), XX(:,2), Membership_values)
% 
% figure()
% scatter3(input_samples(:,1), input_samples(:,2), Membership_values)

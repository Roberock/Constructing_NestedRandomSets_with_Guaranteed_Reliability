function [SVM_save, alphas]=Fuzzy_fit(data)
    
    X = data;
    Y = ones(1,size(X,1));
    [PseudoRadius,PseudoCenter,~] = minimum_boundingSphere(X,1);
    % PseudoCenter=mean(X);
    % PseudoRadius=max(sqrt(sum((X-PseudoCenter).^2,2)).*1.0);
    X_vtmp= randsphere(size(X,1),size(X,2),PseudoRadius.*1.1,'sphere');
    X_virtual=X_vtmp+PseudoCenter;
    Y_virtual = -1*ones(1,size(X_virtual,1));
    Y_fake2class=[Y ,Y_virtual];
    X_fake2class=[X ; X_virtual];

    KernelType= 'rbf';
    
    N_tmp=size(X, 1);
    N = N_tmp;

    i = 1;
    alpha(i) = 0;
    while alpha(i) < 0.95
    
        SVMModel = fitcsvm(X_fake2class,Y_fake2class,'KernelFunction',KernelType,'KernelScale','auto',...
            'BoxConstraint',1,'Verbose',1);
        %SVMModel = fitcsvm(X ,Y ,'KernelFunction',KernelType,'KernelScale','auto',...
        %    'Verbose',1);
        Nsv(i)=sum(SVMModel.IsSupportVector(1:N_tmp));
        
        Beta=10^-4;
        [EpsL(i),EpsU(i)]=epsLU_fast(sum(Nsv(1:i)),N,Beta);
        
        SVM_save{i} = SVMModel;
        alpha(i+1) = EpsU(i);
    
        X_fake2class(SVMModel.IsSupportVector(1:N_tmp),:)=[];
        Y_fake2class(SVMModel.IsSupportVector(1:N_tmp))=[];
        N_tmp=N-Nsv(i);
        i = i + 1;
    end
    
    alphas = EpsU;

end
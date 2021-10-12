function [SupportRemovalStrategy]=Support_removalAlgorithm(MySetFunction,X,Nsets)

Xtmp=X;
N=size(X,1);
Nsv=zeros(Nsets,1);
Eps=zeros(Nsets,2);
SupportIndices=cell(Nsets,1);
for Count=1:Nsets
    [supports,~] = MySetFunction(Xtmp);
    
    if Count==1
        Nsv(Count)=length(supports);
        
    else
        Nsv(Count)=length(supports)+Nsv(Count-1);
        
    end
    SupportIndices{Count}=supports;
    
    Xtmp(supports,:)=[];
    %% remove support constraints
    if  isempty(Xtmp)
        break
    end
    
    %% robustness guarantees
    beta=10^-6;
    [epsL, epsU] = epsLU_fast(Nsv(Count),N,beta); 
    Eps(Count,:)=[epsL, epsU];
    disp(['Nested set number ' num2str(Count) '/' num2str(Nsets) 'computed'])
    %% true reliability
    % use X_test to get true P[x \in Set]
end

SupportRemovalStrategy.Number=Nsv(1:Count);
SupportRemovalStrategy.Indices=SupportIndices;
SupportRemovalStrategy.SetDesignandProp=[];
SupportRemovalStrategy.SetsReliabilityBounds=Eps;
end
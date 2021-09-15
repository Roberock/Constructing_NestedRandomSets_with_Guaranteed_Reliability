N=5000; % number of samples
Nd=2; % we are looking for [min(x), max(x)]
X_fixed=normrnd(1,4,[1,N])';
X_fixed2=X_fixed.^2+normrnd(2,4,[1,N])';
X1=X_fixed;
X2=X_fixed2;
XX=[X1 X2];
%% %  a-posteriori lower and upper bounds
beta=10^-4;
for i=1:N/4
    MaxMin1(i,:)=[min(XX(:,1)), max(XX(:,1))];
    MaxMin2(i,:)=[min(XX(:,2)), max(XX(:,2))];
    XX((XX(:,1)==MaxMin1(i,1)),:)=[];
    XX((XX(:,1)==MaxMin1(i,2)),:)=[];
    XX((XX(:,2)==MaxMin2(i,1)),:)=[];
    XX((XX(:,2)==MaxMin2(i,2)),:)=[];
    [epsL(i), epsU(i)] = epsLU_fast(4*(i-1)+4,N,beta);
end
figure(1)
subplot(2,1,1)
plot(MaxMin1(:,1),epsL,'k')
hold on
plot(MaxMin1(:,2),epsL,'k')
plot(MaxMin1(:,1),epsU,'r')
plot(MaxMin1(:,2),epsU,'r')
subplot(2,1,2)
plot(MaxMin2(:,1),epsL,'k')
hold on
plot(MaxMin2(:,2),epsL,'k')
plot(MaxMin2(:,1),epsU,'r')
plot(MaxMin2(:,2),epsU,'r')
figure (2)
for i=1:50:N/4
    p1=[MaxMin1(i,1) MaxMin2(i,1)];
    p2=[MaxMin1(i,1) MaxMin2(i,2)];
    p3=[MaxMin1(i,2) MaxMin2(i,2)];
    p4=[MaxMin1(i,2) MaxMin2(i,1)];
    Xfillx=[p1(1) p2(1) p3(1) p4(1)];
    Xfilly=[p1(2) p2(2) p3(2) p4(2)];
    h=fill(Xfillx,Xfilly,'r');
    hold on
    set(h,'facealpha',.01)
end?
figure
scatter(MaxMin1(:,1),MaxMin2(:,1),'k'); hold on
scatter(MaxMin1(:,2),MaxMin2(:,2),'k')
scatter(MaxMin1(:,2),MaxMin2(:,1),'k')
scatter(MaxMin1(:,1),MaxMin2(:,2),'k')?
% %% % a-priori upper bound
% X=X_fixed;
% for i=1:N/10
% MaxMin(i,:)=[min(X), max(X)];
%  epsU_apriori(i)  = getepsilon_ConvexDiscard(N,beta,2*(i-1),Nd);
% X(find(X==max(X)))=[];
% X(find(X==min(X)))=[];
% end?
% X_test=normrnd(1,4,[1,10^5])';
% for i=1:2500
% Rel(i)=mean(MaxMin(i,1)< X_test & MaxMin(i,2)> X_test);
% Pf(i)=1-Rel(i);
% end
% plot(MaxMin(:,1),Pf,'b')
% hold on
% plot(MaxMin(:,2),Pf,'b')
% plot(MaxMin(:,1),Rel,':b')
% plot(MaxMin(:,2),Rel,':b')
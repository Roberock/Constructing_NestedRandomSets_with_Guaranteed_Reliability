function fun_plotEclosingSet(X,supports,f_x,varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if nargin <4
    margin_x=0.1;
    margin_y=0.1;
elseif nargin==4
    margin_x=varargin{1};
    margin_y=margin_x;
elseif nargin==5
    margin_x = varargin{1};
    margin_y = varargin{2};
end

%%
[s1,s2]=size(X);
assert(s2==2,'Plot for two dimensions only')

min_X=min(X)-abs(min(X)).*[margin_x,margin_y];
max_X=max(X)+abs(max(X)).*[margin_x,margin_y];
%%
minmaxX = linspace(min_X(1),max_X(1));
minmaxY = linspace(min_X(2),max_X(2));
[MX,MY] = meshgrid(minmaxX,minmaxY);

Z = zeros(100);
for i =1:100
    for j = 1:100
       Z(i,j) = f_x([MX(i,j),MY(i,j)]);
    end
end

%%
figure
contour(MX,MY,Z,[f_x(X(supports(1),:)),f_x(X(supports(1),:))])
hold on 
hold on
scatter(X(:,1),X(:,2)), 
hold on, 
for i=1:length(supports)
scatter(X(supports(i),1),X(supports(i),2),'r')
end
axis equal

end


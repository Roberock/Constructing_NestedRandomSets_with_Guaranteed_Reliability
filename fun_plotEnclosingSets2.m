function fun_plotEnclosingSets2(X,f_x,i_x,varargin)

%% Make some gas particles
N = 5000;
u1=-3+ 6*rand(N,1);
u2=-3+ 20*rand(N,1);

inside=false(N,1); for i =1:N, inside(i)=i_x([u1(i),u2(i)]);end % boolean vactor of particles that are inside the sets

%%
figure, 
scatter(u1(inside),u2inside),'.'), 
hold on
scatter(u1(~inside),u2(~inside),'.')
%%
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

[s1,s2]=size(X);
assert(s2==2,'Plot for two dimensions only')

min_X=min(X)-abs(min(X)).*[margin_x,margin_y];
max_X=max(X)+abs(max(X)).*[margin_x,margin_y];

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
contour(MX,MY,Z,[f_x(X(supports(1),:)),f_x(X(supports(1),:))],'linewidth',2)

scatter(X(:,1),X(:,2),'filled')

for i=1:length(supports)
scatter(X(supports(i),1),X(supports(i),2),'r')
end

axis equal

end

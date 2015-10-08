function X = ACHRsampling(S,b,x,n)
if nargin <4
    n = 100;
else
    n = n+1;
end
nvar = size(S,2);
if nvar ~= length(x)
    fprintf('Column mismatch\n');
    X = x;
    return
end
if size(S,1) ~= length(b)
    fprintf('Constriant mismatch\n');
    X = zeros(nvar,1);
    return
end

%n - number of points
runup = 1000;
discard = 1000;

%initialize empty sample matric
X = zeros(nvar,runup+discard+n);
X(:,1) = x; %initial point

mu = zeros(nvar,1);
sigma = zeros(nvar,1);

%initial run
%generate random direction    
u = rand(nvar,1);
u = u/norm(u); %random direction

%check for intersection between polytope and x0 + ut
Su = S*u;
t = (b-S*x)./Su;

%minimum starting point for line segment within polytope
tmin = max(t(Su<0));

%maximum ending point for line segment within polytope
tmax = min(t(Su>0));

%random point on line segment
yl = tmin-(tmax-tmin)*rand;

%new point
y = x + yl*u;
X(:,2) = y;

%new mean
deltaMu = x-mu;
mu = mu+deltaMu/n;

for i=3:runup+discard+n
    
    %centre to new mean
    %choose a sampled point in random
    v = X(:,randi(i));
    
    %set direction from this point towards mean
    u = (v-mu)/norm(v-mu);
    
    %repeat intersection detection
    Su = S*u;
    t = (b-S*x)./Su;
    tmin = max(t(Su<0));
    tmax = min(t(Su>0));
    X(:,i) = X(:,i-1) + (tmin-(tmax-tmin)*rand)*u;
    
    %new mean 
    deltaMu = X(:,i)-mu;
    mu = mu+deltaMu/n;
end
X = X(:,discard+runup+2:n+discard+runup);
    
% numerical calculation of separatrix from Seydel, 1994
% example using Duffing equation
function w = separatrixSeydel()
% ys1 = [0;0];
% ys2 = [1;0];
% ys3 = [-1;0];
ys1 = [1/4;3/4];
% calculation of jacobian
% J = jacobian(ys1);
J = randomJac(ys1);

% calculation of eigen values
[w,lambda] = eig(J);
lambda = diag(lambda);
w(:,1) = w(:,1)./w(1,1);
w(:,2) = w(:,2)./w(2,2);

% initial points
eps = 1e-4;
z = ys1+eps*w(:,2);

opts = odeset('RelTol',1e-10,'AbsTol',1e-8);
% integrate forward in time
[tf,xaf] = ode45(@randomeq,[0,100],z,opts);
plot(xaf(:,1),xaf(:,2),'r');
hold on

% integrated backward in time
[tb,xab] = ode45(@randomeq,[0,-3],z,opts);
plot(xab(:,1),xab(:,2),'r');

z = ys1-eps*w(:,1);
% integrate forward in time
[tf,xaf] = ode45(@randomeq,[0,100],z,opts);
plot(xaf(:,1),xaf(:,2),'k');
hold on

% integrated backward in time
[tb,xab] = ode45(@randomeq,[0,-40],z,opts);
plot(xab(:,1),xab(:,2),'k');
plot(ys1(1),ys1(2),'LineStyle','none','Marker','.','MarkerSize',20)


end

function dy = duffing(t,y)
% d2y + dy - y + y^3 = 0
% dy1 = y2
% dy2 = y1-y1^3-y2

dy = zeros(2,1);
dy(1) = y(2);
dy(2) = y(1)-y(1)^3-y(2);
end

function dy = randomeq(t,y)
dy = zeros(2,1);
dy(1) = (1-y(1)-y(2))*y(1);
dy(2) = (4-7*y(1)-3*y(2))*y(2);
end

function Jy = jacobian(y)
Jy = [0 1;1-3*y(1)^2 -1];
end

function Jy = randomJac(y)
Jy = [1-2*y(1)-y(2) -y(1);-7*y(2) 4-7*y(1)-6*y(2)];
end

function w = eigenvector(Jy,lambda)
A = Jy-lambda*eye(2);
b = zeros(2,1);
w = A\b;
end
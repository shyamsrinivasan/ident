% [X, Y] = meshgrid(0:0.1:3.5, 0:0.1:1);
% U = X.*(1.5 - 0.5*X - Y);
% V = Y.*(0.75 - 0.125*X - Y);
% L = sqrt((U/3.5).^2 + (V).^2);
% quiver(X, Y, U./L, V./L, 0.5);
% axis tight
% xlabel x
% ylabel y

% [X,Y] = meshgrid(0:0.1:2.5,0:0.1:2.5);
% U = 2.*X.*(1-X./2.5)-1.5.*X.*Y;
% V = -Y+0.8.*X.*Y;
% L = sqrt((U/2.5).^2 + (V./2.5).^2);
% quiver(X, Y, U./L, V./L, 0.5);
% axis tight
% xlabel x
% ylabel y


figure;
for a= 0:0.05:0.5
for b= 0:.1:2
% [t,xa] = ode45(@f, [0,10], [a*b a*(.9/3)*(1-b)]);
% plot(xa(:,1), xa(:,2),'k')
% [t,xa] = ode45(@f, [0,-10], [a*b a*(.9/3)*(1-b)]);
% [t,xa] = ode45(@f, [0,-10], [a;b]);
% plot(xa(:,1), xa(:,2),'r')
end
end
% opts = odeset('RelTol',1e-14,'AbsTol',1e-12);
% [t,xa] = ode45(@fiGEM, [0,1000], [1e6;1e6;1e-6;1e-6],opts);
% plot(t,xa(:,[2 4]));
% axis([0 3.5 0 1])
% hold on
% quiver(X, Y, U./L, V./L, 0.4)
% axis([0 3.5 0 1])


% iGEM example for nullclines
% T = 0.1;
% M1 = 1;
% allx = [];
% for cfp = 0:1:10
%     fun = @(x)0.416.*M1./(1+(cfp./1e4).^2)-(2.5e-4+T).*x;
%     [x,fval,exitflag] = fzero(fun,0.1);
%     allx = [allx;x];
% end

allx = [];
allz = [];
allx3 = [];
a = 0.8;
c = 1;
for y = -2:0.1:2
    fun = @(x)a*x-y;
    [x,fval,exitflag] = fzero(fun,0.01);
    allx = [allx;x];
    fun1 = @(x)y^3-c*x;
    [z,fval,exitflag] = fzero(fun1,0);
    allz = [allz;z];
    fun3 = @(x)x-y;
    [x3,fval,exitflag] = fzero(fun3,0.01);
    allx3 = [allx3;x3];
end

f1 = @(y)y./a;
f2 = @(z)z;
f3 = @(y)y.^3/c;
y = -2:0.1:2;
z = -2:0.1:2;

plot3(f1(y),y,z,'k');
hold on 
plot3(f2(z),y,z,'r');

diff = f1(y)-f2(y);
x1 = f1(y);
x2 = f2(y);
xeq = x1(diff==0);
x = repmat(-1,1,length(y));
plot(y







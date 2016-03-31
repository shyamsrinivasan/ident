M = zeros(4,1);
M(1)  = 1;      % E
M(2)  = 0.01;   % PEP
M(3)  = 0.03;   % FBP
M(4) = 2;       % a.u acetate

dMdt = Kotte_glycolysis(0,M);

% Xss = J*KEXPEP/(vEXmax-J);

[tout,yout] = ode45(@Kotte_glycolysis,0:0.1:200,M);
subplot(221);
plot(tout,yout(:,2));
ylabel('PEP');
subplot(222);
plot(tout,yout(:,3));
ylabel('FBP');
subplot(223);
plot(tout,yout(:,4));
ylabel('acetate');
xlabel('time');
subplot(224);
plot(tout,yout(:,1));
ylabel('super Enzyme E');
xlabel('time');

%phase planes
figure
subplot(131);
plot(yout(:,1),yout(:,2));
xlabel('super Enzyme E');
ylabel('PEP');
subplot(132);
plot(yout(:,1),yout(:,3));
xlabel('super Enzyme E');
ylabel('FBP');
subplot(133);
plot(yout(:,2),yout(:,3));
xlabel('PEP');
ylabel('FBP');


% flux calculation
flux = zeros(length(tout),4);
for it = 1:length(tout)
    flux(it,:) = Kotte_glycolysisflux(yout(it,:));
end
figure
subplot(221);
plot(tout,flux(:,2));
ylabel('Fbp');
subplot(222);
plot(tout,flux(:,3));
ylabel('FBP');
subplot(223);
plot(tout,flux(:,4));
ylabel('E(FBP)');
xlabel('time');
subplot(224);
plot(tout,flux(:,1));
ylabel('super Enzyme E,J');
xlabel('time');

%Bifurcation analysis using MATCONT
% [x,v,s,h,f] = cont(@);


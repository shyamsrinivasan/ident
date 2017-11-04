% script to test Chian et la., 1988's algorithm for determiniing
% approximate regions of attraction
% stable eq (0,0) and saddle eq (1,2)
pt1 = [0;0];
pt2 = [1;2];
fval1 = examplefun(pt1);
fval2 = examplefun(pt2);

fjac1 = examplejac(pt1);
[w1,eig1] = eig(fjac1);
eig1 = diag(eig1);

fjac2 = examplejac(pt2);
[w2,eig2] = eig(fjac2);
eig2 = diag(eig2);

% 2. unstable eigen vectors
unsw2 = w2(:,real(eig2)>0);


% 3. intersection ppoints with e-ball
eps = 1e-3;
% ip_inter = [pt2+eps*unsw2 pt2-eps*unsw2];

figure
hold on
for iw = 1:size(w2,2)
    ip_inter = [pt2+eps*(w2(:,iw)) pt2-eps*(w2(:,iw))];
    itermax = 100;
    saveeps = zeros(size(ip_inter,2),1);
    saveinterpt = zeros(2,2);    
    for ip = 1:size(ip_inter,2)
        % 4. reverse integration of vector field from intersection points
        [~,xrev] = ode45(@odeexamplefun,[0,-1],ip_inter(:,ip));
        iter = 1;
        alpha = 0.9;
        while any(any(abs(xrev-repmat(pt2',size(xrev,1),1))>eps)) && iter<=itermax
    %         alpha = alpha/5;
            eps = eps*alpha;
            if ip==1
                ip_inter(:,ip) = pt2+eps*w2(:,iw);
            elseif ip==2
                ip_inter(:,ip) = pt2-eps*w2(:,iw);
            end 
            [~,xrev] = ode45(@odeexamplefun,[0,-1],ip_inter(:,ip));
            iter = iter+1;
        end
        saveeps(ip) = eps;
        saveinterpt(:,ip) = ip_inter(:,ip);
        % 5. forward intersection from intersection point 
        [tout,xfwd] = ode45(@odeexamplefun,[0,10],saveinterpt(:,ip));
        [tout,xus] =  ode45(@odeexamplefun,[0,-15],saveinterpt(:,ip));
        plot(ip_inter(1,ip),ip_inter(2,ip),'Marker','.','MarkerSize',10,'Color','b');
        plot(xfwd(:,1),xfwd(:,2),'k');
        plot(xus(:,1),xus(:,2),'r');
    end
end




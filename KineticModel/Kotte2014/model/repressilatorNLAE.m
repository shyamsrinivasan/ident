function dX = repressilatorNLAE(x,model,pvec)
% pi = pvec(1:2);
% ps = pvec(3:4);
dX = zeros(length(x),1);
dX = cons(dX,x);

alpha = pvec(1);
beta = pvec(2);
delta = pvec(3);
h = pvec(4);
n = 6;

dX(1) = beta*(x(2)-x(1)); % x0
dX(2) = alpha*((1+delta)/(1+x(5)^h)+delta)-x(2); % y0
dX(3) = beta*(x(4)-x(3)); % x1
dX(4) = alpha*((1+delta)/(1+x(1)^h)+delta)-x(4); % y1
dX(5) = beta*(x(6)-x(5)); % x2
dX(6) = alpha*((1+delta)/(1+x(3)^h)+delta)-x(6); % y2


%  for i = 1:2:n
%      dX(i) = beta*(x(i+1)-x(i));
%      if i>1
%          dX(i+1) = alpha*((1-delta)/(1+x(i-2)^h)+delta)-x(i+1);
%      else
%          dX(i+1) = alpha*((1-delta)/(1+x(n)^h)+delta)-x(i+1);
%      end
%  end


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
for i = 1:2:n
    dX(i) = beta*(x(i+1)-x(i));
    if i>1
        dX(i+1) = alpha*((1-delta)/(1+x(i-2)^h)+delta)-x(i+1);
    else
        dX(i+1) = alpha*((1-delta)/(1+x(n)^h)+delta)-x(i+1);
    end
end


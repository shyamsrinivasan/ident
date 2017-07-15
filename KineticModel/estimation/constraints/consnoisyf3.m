function fx = consnoisyf3(x,p,data)

ratio = 1+x(2)./x(4);
fx = x(6).*(ratio-1).*(ratio).^3 - x(end).*(ratio.^4+p(3).*(1+x(1)./x(5)).^(-4));
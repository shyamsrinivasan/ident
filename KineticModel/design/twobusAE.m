function fx = twobusAE(x,p)
fx = zeros(2,1);
fx(1) = -4*x(2)*sin(x(1)) - p(1);
fx(2) = -4*x(2)^2 + 4*x(2)*cos(x(1)) - p(2);
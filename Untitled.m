o2c = 10; %mM
o2e = 0:0.01:1000; %mM
vo2t = (469.*o2e)./(1+o2e/ko2e+o2c/ko2c);

hold on
plot(o2e,vo2t,'r-');

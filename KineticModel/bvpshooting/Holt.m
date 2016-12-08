% example for shooting with continuation for 2 point BVP
function dy = Holt(t,y)
s = .2;
n = -.1;

dy = zeros(5,1);
dy = cons(dy,y);

dy(1) = y(2);
dy(2) = y(3);
dy(3) = -(3-n)/2*y(1)*y(3) - n*y(2)^2 + 1 - y(4)^2 + s*y(2);
dy(4) = y(5);
dy(5) = -(3-n)/2*y(1)*y(5) - (n-1)*y(2)*y(4) + s*(y(4)-1);
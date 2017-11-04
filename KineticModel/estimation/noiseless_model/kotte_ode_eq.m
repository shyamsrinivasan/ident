function dx = kotte_ode_eq(x,p,flux)

% dx = cell(3,1);

dx = [flux(1)-flux(4)-flux(5);...
      flux(4)-flux(3);...
      flux(2)-flux(6)];

% dx{1} = flux{1}-flux{4}-flux{5};
% dx{2} = flux{4}-flux{3};
% dx{3} = flux{2}-flux{6};

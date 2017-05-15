function fixed_points = kotte_branches(npoints,range,contdir,eqpts,model,pvec,ap,opts)

% initialize parameter
delpar = (range(2) - range(1))/(npoints-1);

% direction of continuation
if contdir == 1
    pararray = range(1):delpar:range(2);
elseif contdir == -1
    pararray = range(2):-delpar:range(1);
end

% continuation of fixed points
[contpt,conteig,contpttype] = continuation(model,pvec,pararray,ap,eqpts,opts);

fixed_points = contpt;
for i = 1:npoints
%     fixed_points = 
end
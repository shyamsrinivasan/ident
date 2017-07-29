% optimize with random initial values
function optsol = multistart_search(prob,x0,optimdata,solveropt)
if nargin<4
%     solveropt = struct('solver','ipopt','multi',1,'multi_pts',[2 2]);
    solveropt = struct('solver','ipopt','multi',0);
end

npts = size(x0,2);
cell_optsol = cell(npts,1);
ts = tic;
parfor ix0 = 1:npts % change this between serial/parallel 
    cell_optsol{ix0} = nlconsopt(prob,x0(:,ix0),solveropt,optimdata);     
end
fprintf('Total time for completion :%4.2g\n',toc(ts));

% convert cell_optsol to structure array
optsol = struct();
for jx0 = 1:npts
    optsol(jx0).x0 = cell_optsol{jx0}.x0;
    optsol(jx0).xval = cell_optsol{jx0}.xval;
    optsol(jx0).fval = cell_optsol{jx0}.fval;
    optsol(jx0).exitflag = cell_optsol{jx0}.exitflag;
    optsol(jx0).info = cell_optsol{jx0}.info;
    optsol(jx0).time = cell_optsol{jx0}.time;
end




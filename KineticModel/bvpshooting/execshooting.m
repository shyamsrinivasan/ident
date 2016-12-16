function [yi,yf,tf,delyi,delyf,flag] =...
execshooting(fh,yi,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,eps,opts,niter)
% execute shooting iteration
if nargin<13
    niter = 100;
end
if nargin<12
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end
if nargin<11
    eps = 1e-4;
end
if nargin<10||isempty(yf)
    yf = zeros(size(yi,1),1);
end
if nargin<9||isempty(delyf)
    delyf = zeros(size(yi,1),1);
end

printBVPstats();
iter = 1;
flag = 1;
while any(abs(delyf(yfknwn))>repmat(eps,length(yfknwn),1))
    tstart = tic;
    if iter<niter
        [yi,yfold,delyi,delyf] = itershooting(fh,yi,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,opts);
        % integrate system with guessed/new intial conditions
        % fprintf('Integrated system to find final value...\n');
        [tfnew,yf,status] = integrateshooting(fh,ti,tf,yi,opts);
        
        if status>0 % successful integration to new final value at new final time
            tf = tfnew;
            % proceed to next iteration
        elseif status<0 % integration fails before new final time at tf
            flag = 0;
            delyf = getvaldiff(yfold,yf);
            printBVPstats(iter,delyf(yfknwn),ti,tf,toc(tstart),'F');   
            tf = [tf tfnew];
            return % restart the BVP solution process till time tf
        end
            
        % check difference in final values
        delyf = getvaldiff(yfold,yf); % assuming ysimf is obtained at tf        
        printBVPstats(iter,delyf(yfknwn),ti,tf,toc(tstart),'S');        
        iter = iter+1;
    else
        fprintf('Number of iterations exceeded\n');
        flag = 0;
    end
end



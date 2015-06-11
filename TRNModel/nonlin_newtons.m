%Newtons method to determine roots of Nonlinear vector valued function
function [xold] = nonlin_newtons(initval,defparval,ng,model,maxiter)
% Formulate different function for this Newtons using cons() from ADMAT
call_func = @(X)admat_intODEmodel(X,defparval,ng,model);
% initial evaluation
xold = initval;
V = eye(sum(ng)-ng(4));
iter = 1;
fprintf('Iteration  Error   Function Error');
while iter <= maxiter
    admat_obj = deriv(xold,V);
    admat_result = call_func(admat_obj);
    Fold = getval(admat_result);
    Jold = getydot(admat_result);
    xnew = xold - Jold\Fold;   
    Error = norm(xold-xnew);
    Fnew = getval(call_func(deriv(xnew,V)));
    Ferror = norm(Fnew);
    fprintf('%d     %6.4g     %6.4g     \n',iter,Error,Ferror);
    if Error <= 1e-6
        fprintf('\nDesired tolerance achieved\nStopping Iteration at %d\n',iter);
        iter = maxiter;        
    else
        xold = xnew;
        iter = iter + 1;
    end
end
xold = xnew;

        
%x(t+1) = x(t) - J(x(t))\F(x(t))
%stop condition = norm(f(x(t+1))) < tol1
%or
%stop condition = max(abs(x(t)-x(t+1)))/max(abs(x(t+1))) < tol2

% nexpt = length(fieldnames(allSolution));
% x = zeros(nexpt,1);
% y = zeros(nexpt,1);
% for iexpt = 1:nexpt
%     sim_name = sprintf('sim_%d',iexpt);
%     x(iexpt) = allfinalSS.(sim_name).y(1);  
%     y(iexpt) = allfinalSS.(sim_name).y(13);    
% end
% plot(x,y);
% xlabel('mRNA umole/g DCW');
% ylabel('Protein umole/g DCW');
return

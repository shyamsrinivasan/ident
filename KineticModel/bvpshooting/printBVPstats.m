function printBVPstats(iter,delyf,t0,tf)
if nargin<1
    call=1;
else
    call=2;
end

if call == 1
    fprintf('\n\tTwo point BVP solution stats\t\t\n');
    fprintf('--------------------------------------\n');
    fprintf('iter \t\t max(delyf) \t\t t0 \t\t tf \n');
elseif call == 2
    fprintf('%d \t\t\t %5.3g \t\t\t %3.2f \t\t\t %3.2f \n',iter,max(delyf),t0,tf);
end
% collect data s.t. x-axis is experimental data
% yaxis is opt and calc data
function plot_data = scatterdata(nplots,exp_xss,exp_xerr,opt_xss,opt_xerr,calc_xss,calc_xerr)
if nargin<7
    calc_xerr = [];    
end
if nargin<6
    calc_xss = [];    
end
calc = 1;opt = 1;
if isempty(calc_xss)
   calc = 0;
end
if isempty(opt_xss)
   opt = 0;
end
if calc && isempty(calc_xerr)
    calc_xerr = zeros(size(calc_xss,1),size(calc_xss,2));
end
if opt && isempty(opt_xerr)
    opt_xerr = zeros(size(opt_xss,1),size(opt_xss,2));
end
if isempty(exp_xerr)
    exp_xerr = zeros(size(exp_xss,1),size(exp_xss,2));
end


plot_data.x
plot_data.y





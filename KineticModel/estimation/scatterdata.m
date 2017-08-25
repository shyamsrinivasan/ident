% collect data s.t. x-axis is experimental data
% yaxis is opt and calc data
function plot_data =...
scatterdata(nplots,exp_xss,exp_xerr,opt_xss,opt_xerr,calc_xss,calc_xerr,type)
if nargin<8
    type=1; % averages
end
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

if type==1
    nvar = size(exp_xss,1);
    ss_x = cell(nvar,1);
    ss_y = cell(nvar,2);
    err_x = cell(nvar,1);
    err_y = cell(nvar,2);
    for j = 1:nplots
        ss_x{j} = exp_xss(j,:)';
        err_x{j} = exp_xerr(j,:)';
        if opt && calc        
            ss_y{j,1} = opt_xss(j,:)';
            err_y{j,1} = opt_xerr(j,:)';
            ss_y{j,2} = calc_xss(j,:)';        
            err_y{j,2} = calc_xerr(j,:)';        
        elseif opt && ~calc
            ss_y{j,2} = [];
            err_y{j,2} = [];
            ss_y{j,1} = opt_xss(j,:)';
            err_y{j,1} = opt_xerr(j,:)';
        elseif ~opt && calc
            ss_y{j,1} = [];
            err_y{j,1} = [];
            ss_y{j,2} = calc_xss(j,:)';        
            err_y{j,2} = calc_xerr(j,:)';
        end
    end
elseif type==2 % parse through data
    
end

plot_data.x = ss_x;
plot_data.y = ss_y;
plot_data.x_err = err_x;
plot_data.y_err = err_y;





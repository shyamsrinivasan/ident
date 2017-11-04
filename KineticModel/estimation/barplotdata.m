function plot_data =...
barplotdata(nplots,exp_xss,exp_xerr,opt_xss,opt_xerr,calc_xss,calc_xerr)
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

nvar = size(exp_xss,1);
% use cell arrays to store data for each nvar subplot
ss_plot = cell(nvar,1);
ss_error = cell(nvar,1);

% collect all data
for j = 1:nplots
    if opt && calc
        ss_plot{j} = [exp_xss(j,:)' opt_xss(j,:)' calc_xss(j,:)'];
        ss_error{j} = [exp_xerr(j,:)' opt_xerr(j,:)' calc_xerr(j,:)'];
    elseif opt && ~calc
        ss_plot{j} = [exp_xss(j,:)' opt_xss(j,:)'];
        ss_error{j} = [exp_xerr(j,:)' opt_xerr(j,:)'];
    elseif ~opt && calc
        ss_plot{j} = [exp_xss(j,:)' calc_xss(j,:)'];
        ss_error{j} = [exp_xerr(j,:)' calc_xerr(j,:)'];
    end    
end

plot_data.ss = ss_plot;
plot_data.err = ss_error;

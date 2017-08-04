% function to plot 2-d and 3-d diagrams comapring parameters
function compare_pars(est_data)
if isempty(est_data)
    fprintf('No optimal solution found\n');
    return
end

figure
par = est_data.xpar;
[npar,~] = size(par);

if npar<=2 % 2-d plot
    line(par(1,:),par(2,:),'LineStyle','none',...
                            'Marker','.',...
                            'MarkerSize',10);
elseif npar<=3 % 3-d plot
    line(par(1,:),par(2,:),par(3,:),'LineStyle','none',...
                                    'Marker','.',...
                                    'MarkerSize',10);
end




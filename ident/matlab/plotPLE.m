% function to plot profile likelihood curves with thresholds
% ci_df1 - point wise CI
% ci_dfall - full dof CI
function plotPLE(PLEval,ci_df1,ci_df_all)
if nargin<3
    ci_df_all = [];
end
if nargin<2
    ci_df1 = [];
end

if isfield(PLEval,'thetai_inc')
    thetai = PLEval.thetai_inc;
end
if isfield(PLEval,'chiPLE')
    chi2PLE = PLEval.chiPLE;
end

npts = length(thetai);
figure
line(log10(thetai),chi2PLE,'LineStyle','-','LineWidth',1.5,'Color','k');
hold on
if ~isempty(ci_df1)
    line(log10(thetai),repmat(ci_df1,1,npts),'LineStyle','-.',...
                                            'LineWidth',2,...
                                            'Color','k');
end
if ~isempty(ci_df_all)
    line(log10(thetai),repmat(ci_df_all,1,npts),'LineStyle','-.',...
                                                'LineWidth',2,...
                                                'Color','k');
end

% function to plot profile likelihood curves with thresholds
% ci_df1 - point wise CI
% ci_dfall - full dof CI
function plotPLE(PLEval,ci_df1,ci_df_all)

if isfield(PLEval,'thetai_inc')
    thetai = PLEval.thetai_inc;
end
if isfield(PLEeval,'chi2PLE')
    chi2PLE = PLEval.chi2PLE;
end

npts = length(thetai);
figure
line(log10(thetai),chi2PLE,'LineStyle','-','LineWidth',1.5,'Color','k');
hold on
line(log10(thetai),repmat(ci_df1,1,npts),'LineStyle','-.','LineWidth',2,'Color','k');
line(log10(thetai),repmat(ci_df_all,1,npts),'LineStyle','-.','LineWidth',2,'Color','k');

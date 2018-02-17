function plot2d_dynamic(solution, concentration, flux)
if nargin<3
    flux=0;
end
if nargin<2
    concentration=0;
end
nplots = size(solution,2);
if nplots>1    
    ncolumns = 2;
    if rem(nplots,2)==0
        nrows = nplots/2;
    else
        nrows = (nplots+1)/2;
    end    
else
    ncolumns = 1;
    nrows = 1;
end
haxis = zeros(nplots, 1);

if concentration
    hfig = figure;
    plot_concentration(solution,hfig,haxis,nplots,nrows,ncolumns);
end
if flux
    hfig = figure;
    plot_flux(solution,hfig,haxis,nplots,nrows,ncolumns);
end

    
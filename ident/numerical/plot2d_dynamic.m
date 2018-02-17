function plot2d_dynamic(solution, concentration, flux)
if nargin<3
    flux=0;
end
if nargin<2
    concentration=0;
end
hfig = figure;
nplots = size(solution,2);
ncolumns = 2;
if nplots>1    
    if rem(ndata,2)==0
        nrows = ndata/2;
    else
        nrows = (ndata+1)/2;
    end    
else
    nrows = 1;
end
haxis = zeros(nplots, 1);

if concentration
    plot_concentration(solution,hfig,haxis,nplots,nrows,ncolumns);
end
if flux
    plot_flux(solution,hfig,haxis,nplots,nrows,ncolumns);
end

    
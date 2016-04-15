function new_data = MonitorCVodeGraph(call,T,Y,YQ,YS,data)
if nargin ~= 6
    error('Mointor data undefined');
end
new_data = [];
if call == 0
    %initialize uninitialzed values
    data = initialize(data);
    if data.post
        if data.grph
            if data.stats | data.cntr
                data.hfg = figure;
            end
            %     Number of subplots in figure hfg
            if data.stats
                data.npg = data.npg + 2;
            end
            if data.cntr
                data.npg = data.npg + 1;
            end
        end
    end
    
    %initialize other data
    data.i = 0;
    data.n = 1;
    data.t = zeros(1,data.updt);
    if data.stats
        data.h = zeros(1,data.updt);
        data.q = zeros(1,data.updt);
    end
    if data.cntr
        data.nst = zeros(1,data.updt);
        data.nfe = zeros(1,data.updt);
        data.nni = zeros(1,data.updt);
        data.netf = zeros(1,data.updt);
        data.ncfn = zeros(1,data.updt);
    end
    data.first = true;
    data.initialized = false;
    new_data = data;
    return
else
    if data.first
    end
    % Extract variables from data
    hfg = data.hfg;
    npg = data.npg;
    i = data.i;
    n = data.n;
    t = data.t;
    h = data.h;
    q = data.q;
    nst = data.nst;
    nfe = data.nfe;
    nni = data.nni;
    netf = data.netf;
    ncfn = data.ncfn;
end
%Load current stats
if call == 1
    if i~=0
        i = i-1;
        data.i = i;
        new_data = data;
        return
    end
    si = CVodeGetStats;
    t(n) = si.tcur;
    if data.stats
        h(n) = si.hlast;
        q(n) = si.qlast;
    end
    if data.cntr
        nst(n) = si.nst;
        nfe(n) = si.nfe;
        nni(n) = si.nni;
        netf(n) = si.netf;
        ncfn(n) = si.ncfn;
    end
end
%check if time to post
if data.post & (n==data.updt|call==2)
    if call==2
        n=n-1;
    end
    if ~data.initialized
        %initialize data
        graph_init(n,hfg,npg,data.stats,data.cntr,t,h,q,nst,nfe,nni,netf,ncfn);
        data.initialized = true;        
    else
        %update data
        graph_update(n,hfg,npg,data.stats,data.cntr,t,h,q,nst,nfe,nni,netf,ncfn);
    end
    
    if call==2
        %final data output
        graph_final(hfg,npg,data.cntr,data.stats);
        return
    end
    n = 1;
else
    n = n+1;
end
%save updated values in data
data.i = data.skip;
data.n = n;
data.npg = npg;
data.h = h;
data.q = q;
data.nst = nst;
data.nfe = nfe;
data.nni = nni;
data.netf = netf;
data.ncfn = ncfn;
new_data = data;
return

function data = initialize(data)
if ~isfield(data,'mode')
    data.mode = 'graphical';
end
if ~isfield(data,'updt')
    data.updt = 0;
end
if ~isfield(data,'skip')
    data.skip = 0;
end
if ~isfield(data,'stats')
  data.stats = true;
end
if ~isfield(data,'cntr')
  data.cntr = true;
end
if ~isfield(data,'post')
  data.post = true;
end
data.grph = true;

%other initializations
data.npg = 0;
data.hfg = 0;
data.h = 0;
data.q = 0;
data.nst = 0;
data.nfe = 0;
data.nni = 0;
data.netf = 0;
data.ncfn = 0;

function graph_init(n,hfg,npg,stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)
fig_name = 'ODE run stats';

figure(hfg);
set(hfg,'Name',fig_name);
set(hfg,'color',[1 1 1]);
pl = 0;

tlab = 'time';
%step size and order
if stats
    pl = pl+1;
    subplot(npg,1,pl);
    semilogy(t(1:n),abs(h(1:n)),'-');
    hold on 
    box on
    grid on
    xlabel(tlab);
    ylabel('|step size|');
    
    pl = pl+1;
    subplot(npg,1,pl);
    plot(t(1:n),q(1:n),'-');
    hold on 
    box on 
    grid on
    xlabel(tlab);
    ylabel('Order');
end

%counters
if cntr
    pl = pl+1;
    subplot(npg,1,pl);
    plot(t(1:n),nst(1:n),'k-');
    hold on
    plot(t(1:n),nfe(1:n),'b-');
    plot(t(1:n),nni(1:n),'r-');
    plot(t(1:n),netf(1:n),'g-');
    plot(t(1:n),ncfn(1:n),'c-');
    hold on
    box on
    grid on
    xlabel(tlab);
    ylabel('Counters');
end
drawnow

function graph_update(n,hfg,npg,stats,cntr,t,h,q,nst,nfe,nni,netf,ncfn)
figure(hfg);
pl = 0;

%step size and order
if stats
    pl = pl+1;
    subplot(npg,1,pl);
    hc = get(gca,'Children');
    xd = [get(hc,'Xdata') t(1:n)];
    yd = [get(hc,'YData') abs(h(1:n))];
    set(hc,'XData',xd,'YData',yd);
    
    pl = pl+1;
    subplot(npg,1,pl);
    hc = get(gca,'Children');
    xd = [get(hc,'Xdata') t(1:n)];
    yd = [get(hc,'YData') q(1:n)];
    set(hc,'XData',xd,'YData',yd);
end

%counters
if cntr
    pl = pl+1;
    subplot(npg,1,pl);
    hc = get(gca,'Children');
    % Attention: Children are loaded in reverse order!
    xd = [get(hc(1),'XData') t(1:n)];
    yd = [get(hc(1),'YData') ncfn(1:n)];
    set(hc(1), 'XData', xd, 'YData', yd);
    yd = [get(hc(2),'YData') netf(1:n)];
    set(hc(2), 'XData', xd, 'YData', yd);
    yd = [get(hc(3),'YData') nni(1:n)];
    set(hc(3), 'XData', xd, 'YData', yd);
    yd = [get(hc(4),'YData') nfe(1:n)];
    set(hc(4), 'XData', xd, 'YData', yd);
    yd = [get(hc(5),'YData') nst(1:n)];
    set(hc(5), 'XData', xd, 'YData', yd);
end
drawnow

function graph_final(hfg,npg,stats,cntr)
figure(hfg);
pl = 0;

if stats
    pl = pl+1;
    subplot(npg,1,pl);
    hc = get(gca,'Children');
    xd = get(hc,'XData');    
    if min(xd) ~= max(xd)
        set(gca,'XLim',[min(xd) max(xd)]);
    end
    
    pl = pl+1;
    subplot(npg,1,pl);
    ylim = get(gca,'YLim');
    ylim(1) = ylim(1)-1;
    ylim(2) = ylim(2) + 1;
    set(gca,'YLim',ylim);
    if min(xd)~=max(xd)
        set(gca,'XLim',[min(xd) max(xd)]);
    end
end
if cntr
    pl = pl+1;
    subplot(npg,1,pl);
    hc = get(gca,'Children');
    xd = get(hc(1),'XData');
    if min(xd)~=max(xd)
        set(gca,'XLim',[min(xd) max(xd)]);
    end
    legend('nst','nfe','nni','netf','ncfn',2);
end
    
        
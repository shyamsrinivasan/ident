%New ODE Monitor Function for both TRN and Kinetic Models
function [newdata] = getSolverStats(call,T,Y,YQ,YS,data)
newdata = [];

if call == 0
    %Initialize non-existing fields
    data = initialize_data(data);
    if data.post
        if data.stats
            %Create Figure
            data.hfg = figure;
        end
        %Subplots for each statistic plotted
        if data.stats
            data.npg = data.npg + 2;
        end   
    end
    
    if data.stats
        data.hlast = zeros(1,data.update);
        data.qlast = zeros(1,data.update);        
    end
    if data.cntr
        data.nst = zeros(1,data.update);
        data.nfe = zeros(1,data.update);
        data.nni = zeros(1,data.update);
        data.netf = zeros(1,data.update);
        data.nctf = zeros(1,data.update);
    end
    data.n = 1;
    data.t = zeros(1,data.update);
    data.first = true;
    data.initialze = false;%since post = false => graph windows are uninitialized
    newdata = data;
    return;
        
else %If not initialization call (call ~= 0)
    if data.first %if first call
        %Initialize to print solution(Y), sensitivity(YS)
        data.first = false;
    end        
    %Extract Info from data
    hfg = data.hfg;
    npg = data.npg;
    hlast = data.hlast;    
    qlast = data.qlast;    
    nst = data.nst;
    nfe = data.nfe;
    nni = data.nni;
    netf = data.netf;
    nctf = data.nctf;
    n = data.n;
    t = data.t;
end

if call == 1
    %Load current stats
    stats = CVodeGetStats;
    t(n) = stats.tcur;
    if data.stats
        hlast(n) = stats.hlast;        
        qlast(n) = stats.qlast;        
    end
    if data.cntr
        nst(n) = stats.nst;
        nfe(n) = stats.nfe;
        nni(n) = stats.nni;
        netf(n) = stats.netf;
        nctf(n) = stats.ncfn;
    end
end

if data.post && (n == data.update || call == 2)
    
    if call == 2
        n = n-1;
    end
    if ~data.initialized
        %Initialize with new figure/text
        if data.stats
            graph_init(n,hfg,npg,data.stats,t,hlast,qlast);
        end
        if data.cntr
            print_text(n,data.cntr,t,nst,nfe,nni,netf,nctf);
        end
    else
        %Update existing figures/texts
        if data.stats
            graph_update(n,hfg,npg,data.stats,t,hlast,qlast);
        end
        if data.cntr
            print_text(n,data.cntr,t,nst,nfe,nni,netf,nctf);
        end
    end
    if call == 2
        %Print final values
        if data.stats
            graph_final(hfg,npg,stats);
        end
        if data.cntr
            print_text(n,data.cntr,t,nst,nfe,nni,netf,nctf);
        end
        return;
    end
    n = 1;
else
    n= n + 1;
end
%Update values
data.n = n;
data.npg = npg;
data.t = t;
data.hlast = hlast;
data.qlast = qlast;
data.nst  = nst;
data.nfe  = nfe;
data.nni  = nni;
data.netf = netf;
data.nctf = ncfn;
newdata = data;
return;
    
function data = initialize_data(data)
    %fields
    if ~isfield(data,'update')
        data.update = 50;
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

    %Counters
    data.npg = 0;
    data.hfg = 0;
    data.hlast = 0;
    data.qlast = 0;
    data.nst = 0;
    data.nfe = 0;
    data.nni = 0;
    data.netf = 0;
    data.ncfn = 0;      
end
    
function graph_init(n,hfg,npg,stats,t,hlast,qlast)
    figname = 'Run Statistics';
    figure(hfg);
    set(hfg,'Name',figname);
    set(hfg,'color',[1 1 1]);
    pl = 0;

    tlabel = 'Time \rightarrow';
    if stats
        pl = pl + 1;
        subplot(npg,1,pl);
        semilogy(t(1:n),abs(hlast(1:n)),'-');
        hold on;
        box on;
        grid on;
        xlabel(tlabel);
        ylabel('|Step Size|');

        pl = pl + 1;
        subplot(npg,1,pl);
        plot(t(1:n),qlast(1:n),'-');
        hold on;
        box on;
        grid on;
        xlabel(tlabel);
        ylabel('Method Order');
    end
    drawnow;
end

function graph_update(n,hfg,npg,stats,t,hlast,qlast)
    figure(hfg);
    pl = 0;
    if stats
        pl = pl + 1;
        subplot(npg,1,pl);
        hc = get(gca,'Children');
        newx = [get(hc,'XData') t(1:n)];
        newy = [get(hc,'YData') abs(hlast(1:n))];
        set(hc,'XData',newx,'YData',newy);

        pl = pl + 1;
        subplot(npg,1,pl);
        hc = get(gca,'Children');
        newx = [get(hc,'XData') t(1:n)];
        newy = [get(hc,'YData') qlast(1:n)];
        set(hc,'XData',newx,'YData',newy);
    end
    drawnow;
end

function graph_final(hfg,npg,stats)
    figure(hfg)
    pl = 0;        
    if stats
        pl = pl + 1;
        subplot(npg,1,pl);
        hc = get(gca,'Children');
        newx = get(hc,'XData');
        set(gca,'XLim',sort([newx(1) newx(end)]));

        pl = pl + 1;
        subplot(npg,1,pl);
        ylim = get(gca,'YLim');
        set(gca,'YLim',[ylim(1)-1 ylim(2)+1]);
        set(gca,'XLim',sort([newx(1) newx(end)]));
    end
end
        
function print_text(n,cntr,t,nst,nfe,nni,netf,nctf)
    %Print Run Stats to screen (commnad prompt)
    if cntr
        for i = 1:n
            fprintf('%10.3e   %5d   %5d   %5d   %5d   %5d  \n',...
                     t(i),nst(i),nfe(i),nni(i),netf(i),nctf(i)); 
        end
    end
end      
end
    
    
%Following are printed directly to screen at the final call (call == 2)
%nst - number of integration steps
%nfe - number of right hand side function calls
%netf - number of error test failures
%ncfn - number of convergence test failures
%nni - number of nonlinear solver iterations

%Following are plotted/updated at every function call (call == 1)
%tcur     - current time reached by the integrator
%hlast    - last step size used
%qlast - last method order used


    
    %Initialize
%   It is called after every internal CVode step and can be used to
%   monitor the progress of the solver. MONFUN is called with CALL=0
%   from CVodeInit at which time it should initialize itself and it
%   is called with CALL=2 from CVodeFree. Otherwise, CALL=1.
%
%   It receives as arguments the current time T, solution vector Y,
%   and, if they were computed, quadrature vector YQ, and forward 
%   sensitivity matrix YS. If YQ and/or YS were not computed they
%   are empty here.
%
%   If additional data is needed inside MONFUN, it must be defined
%   as
%      FUNCTION NEW_MONDATA = MONFUN(CALL, T, Y, YQ, YS, MONDATA)
%   If the local modifications to the user data structure need to be 
%   saved (e.g. for future calls to MONFUN), then MONFUN must set
%   NEW_MONDATA. Otherwise, it should set NEW_MONDATA=[] 
%   (do not set NEW_MONDATA = DATA as it would lead to unnecessary copying).
%
%   A sample monitoring function, CVodeMonitor, is provided with CVODES.
%
%   See also CVodeSetOptions, CVodeMonitor
%
%   NOTES: 
%   
%   MONFUN is specified through the MonitorFn property in CVodeSetOptions. 
%   If this property is not set, or if it is empty, MONFUN is not used.
%   MONDATA is specified through the MonitorData property in CVodeSetOptions.
%
%   See CVodeMonitor for an implementation example.

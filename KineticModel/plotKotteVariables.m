function plotKotteVariables(xdata,ydata,cnt)
figure
switch cnt
    case 1 % concentrations
        [~,~,np] = size(ydata);
        if np>1
            for ip = 1:np
                subplot(2,2,1);
                plot(xdata,ydata(:,1,ip));
                ylabel('super Enzyme E');
                hold on
                subplot(2,2,2);
                plot(xdata,ydata(:,2,ip));
                ylabel('PEP');
                hold on
                subplot(2,2,3);
                plot(xdata,ydata(:,3,ip));
                ylabel('FBP');
                xlabel('time');
                hold on
            end
        else
            subplot(2,2,1);
            plot(xdata,ydata(:,1));
            ylabel('super Enzyme E');    
            subplot(2,2,2);
            plot(xdata,ydata(:,2));
            ylabel('PEP');    
            subplot(2,2,3);
            plot(xdata,ydata(:,3));
            ylabel('FBP');
            xlabel('time');    
        end
    case 2 % fluxes
        [~,~,np] = size(ydata);
        if np>1            
            for ip = 1:np
                subplot(2,2,1);
                plot(xdata,ydata(:,1,ip));
                ylabel('super Enzyme E,J');
                hold on
                subplot(2,2,2);
                plot(xdata,ydata(:,2,ip));
                ylabel('vEX');
                hold on
                subplot(2,2,3);
                plot(xdata,ydata(:,3,ip));
                ylabel('vFbp');
                xlabel('time');
                hold on
                subplot(2,2,4)
                plot(xdata,ydata(:,4,ip));
                ylabel('E(FBP)');
                xlabel('time');
                hold on
            end
        else            
            subplot(2,2,1);
            plot(xdata,ydata(:,1));
            ylabel('super Enzyme E,J');
            subplot(2,2,2);
            plot(xdata,ydata(:,2));
            ylabel('vEX');
            subplot(2,2,3);
            plot(xdata,ydata(:,3));
            ylabel('vFbp');
            xlabel('time');  
            subplot(2,2,4)
            plot(xdata,ydata(:,4));
            ylabel('E(FBP)');
            xlabel('time');
        end
    case 3 % scatter 
        subplot(2,2,1)
        plot(xdata(1,:),ydata(1,:),'LineStyle','none','Marker','o','MarkerSize',5);
        ylabel('super Enzyme E');
        subplot(2,2,2);
        plot(xdata(2,:),ydata(2,:),'LineStyle','none','Marker','o','MarkerSize',5);
        ylabel('PEP');
        subplot(2,2,3);
        plot(xdata(3,:),ydata(3,:),'LineStyle','none','Marker','o','MarkerSize',5);
        ylabel('FBP');
end

    
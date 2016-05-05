function plotKotteVariables(xdata,ydata,cnt)
figure
switch cnt
    case 1 % concentrations
        [~,~,np] = size(ydata);
        if np>1
            for ip = 1:np
                subplot(2,2,1);
                plot(xdata,ydata(:,1,ip),'Color','k','LineWidth',2);
                ylabel('super Enzyme E');
                hold on
                subplot(2,2,2);
                plot(xdata,ydata(:,2,ip),'Color','k','LineWidth',2);
                ylabel('PEP');
                hold on
                subplot(2,2,3);
                plot(xdata,ydata(:,3,ip),'Color','k','LineWidth',2);
                ylabel('FBP');
                xlabel('time');
                hold on
            end
        else
            subplot(2,2,1);
            plot(xdata,ydata(:,1),'Color','k','LineWidth',2);
            ylabel('super Enzyme E');    
            subplot(2,2,2);
            plot(xdata,ydata(:,2),'Color','k','LineWidth',2);
            ylabel('PEP');    
            subplot(2,2,3);
            plot(xdata,ydata(:,3),'Color','k','LineWidth',2);
            ylabel('FBP');
            xlabel('time');    
        end
    case 2 % fluxes
        [~,~,np] = size(ydata);
        if np>1            
            for ip = 1:np
                subplot(2,3,1);
                plot(xdata,ydata(:,1,ip),'Color','k','LineWidth',2);
                ylabel('super Enzyme E,J');
                hold on
                subplot(2,3,2);
                plot(xdata,ydata(:,2,ip),'Color','k','LineWidth',2);
                ylabel('E(FBP)');                
                hold on
                subplot(2,3,3);
                plot(xdata,ydata(:,3,ip),'Color','k','LineWidth',2);
                ylabel('vFbp');
                xlabel('time');
                hold on
                subplot(2,3,4)
                plot(xdata,ydata(:,4,ip),'Color','k','LineWidth',2);
                ylabel('vEX');
                xlabel('time');
                hold on
                subplot(2,3,5)
                plot(xdata,ydata(:,5,ip),'Color','k','LineWidth',2);
                ylabel('PEPout');
                xlabel('time');
                hold on
            end
        else            
            subplot(2,3,1);
            plot(xdata,ydata(:,1),'Color','k','LineWidth',2);
            ylabel('super Enzyme E,J');
            subplot(2,3,2);
            plot(xdata,ydata(:,2),'Color','k','LineWidth',2);
            ylabel('E(FBP)');
            subplot(2,3,3);
            plot(xdata,ydata(:,3),'Color','k','LineWidth',2);
            ylabel('vFbp');
            xlabel('time');  
            subplot(2,3,4)
            plot(xdata,ydata(:,4),'Color','k','LineWidth',2);
            ylabel('vEX');            
            xlabel('time');
            subplot(2,3,5);
            plot(xdata,ydata(:,5),'Color','k','LineWidth',2);
            ylabel('PEPout');
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

    
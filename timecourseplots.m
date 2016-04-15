function timecourseplots(time,y,type,model,name)
if nargin<5
    name = {};
end
switch(type)
    case 'c'
        metCurves(name,y)
    case 'f'
        fluxCurves(name,y)
    case 'both'
        metCurves(name,y)
        fluxCurves(name,y)
    otherwise
        fprintf('No suitable plot functions found\n');
end

function metCurves(name,y)
    if isempty(name)
        name = model.mets(1:model.nint_metab);
        fprintf('Printing all Metabolites in a single plot\n');
    else
        name = name.mets;
    end
    nplots = length(name);
    if nplots == model.nint_metab
        hfig = figure('Name','metabolites');   
        figure(hfig);
        plot(gca,time,y(1:model.nint_metab,:));
%         plotData(time,y(1:model.nint_metab,:),hfig);
    else
        hst = zeros(nplots,1);
        for in = 1:nplots
            tfm = strcmpi(model.mets,name{in});
            if any(tfm)
                % create new figure
                if isempty(findobj('type','figure','Name','metabolites'))
                    hfig = figure('Name','metabolites');
                else
                    hfig = findobj('type','figure','Name','metabolites');
                end
                if rem(nplots,2)==0
                    hst(in) = subplot(nplots/2,2,in);
                else
                    hst(in) = subplot((nplots+1)/2,2,in);
                end
                LineP.DisplayName = model.mets{tfm};
                LineP.YLabel = model.mets{tfm};
                plotData(time,y(tfm,:),[hfig,hst(in)],LineP);                
            end
        end
    end
end

function fluxCurves(name,y)
    if isempty(name)
        name = model.rxns(model.Vind);
        fprintf('Printing all fluxes in a single plot\n');
    else
        name = name.rxns;
    end
    nplots = length(name);
    if nplots == length(model.Vind)
        hfig = figure('Name','fluxes');   
        figure(hfig);
        plot(gca,time,y);
%         plotData(time,y(model.Vind,:),[hfig,hst]);
    else
        hst = zeros(nplots,1);
        for in = 1:nplots
            tff = strcmpi(model.rxns,name{in});
            if any(tff)
                % create new figure
                if isempty(findobj('type','figure','Name','fluxes'))
                    hfig = figure('Name','fluxes');
                else
                    hfig = findobj('type','figure','Name','fluxes');
                end
                if rem(nplots,2)==0
                    hst(in) = subplot(nplots/2,2,in);
                else
                    hst(in) = subplot((nplots+1)/2,2,in);
                end
                LineP.DisplayName = model.rxns{tff};
                LineP.YLabel = model.rxns{tff};
                plotData(time,y(tff,:),[hfig,hst(in)],LineP);    
            end
        end
    end
end

end

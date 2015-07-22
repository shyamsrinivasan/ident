function flux = ExFlux(model,Y,flux,Vex,type,EC)
%Calculate exchange fluxes in and outside the ODE
useVmax = 0;
if nargin < 6
    useVmax = 1;
end
k = .1;
K = 1;
Kmm = 100;
if useVmax
    for iv = 1:length(Vex)    
        subs = model.S(:,Vex(iv))<0;
        prud = model.S(:,Vex(iv))>0;
        if ~isempty(subs) && ~isempty(prud)
            if strcmpi(type,'mm')
                flux(Vex(iv)) = k*prod(Y(subs))/(Kmm+prod(Y(prud)));
            else
                if model.rev(Vex(iv))
                    flux(Vex(iv)) = k*prod(Y(subs)) - K*prod(Y(prud));
                else
                    flux(Vex(iv)) = k*prod(Y(subs));
                end      
            end
        end
    end
else
    for iv = 1:length(Vex)    
        subs = model.S(:,Vex(iv))<0;
        prud = model.S(:,Vex(iv))>0;
        if ~isempty(subs) && ~isempty(prud)
            if strcmpi(type,'mm')
                flux(Vex(iv)) = k*prod(Y(subs))/(Kmm+prod(Y(prud)));
            else
                if model.rev(Vex(iv))
                    flux(Vex(iv)) = (k*prod(Y(subs)) - K*prod(Y(prud)));
                    %*EC(Vex(iv));
                else
                    flux(Vex(iv)) = k*prod(Y(subs));%*EC(Vex(iv));
                end 
            end
        end
    end
end
return
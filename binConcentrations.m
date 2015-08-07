%Define the uiqueness of concentrations
function [conc,multiSS] = binConcentrations(conc)
[nvar,nsampl] = size(conc(:,2:end));
% eps = 1e-5;
conc(:,1) = roundsd(conc(:,1),5,'ceil');
multiSS = cell(nvar,1);
% newmultiSS = cell(nvar,1);
for iv = 1:nvar
    nci = conc(iv,2:end);    
    nciint = floor(nci);
    nciflt = roundsd(nci-nciint,3,'ceil');
    nci = nciint+nciflt;
    
    for ism = 1:nsampl        
        if ism < nsampl
            if nci(1,ism+1) > 1
                if abs(nci(1,ism)-nci(1,ism+1))<= 1e0
                    nci(1,ism+1) = nci(1,ism);              
                end
            else
                if abs(nci(1,ism)-nci(1,ism+1))<= 1e-6
                    nci(1,ism+1) = nci(1,ism);              
                end
            end
        end
    end    
    
    conc(iv,2:end) = nci;
    tab = tabulate(nci);    
    ncnt = sum(tab(:,2));
    iva = 1;    
    notiva = tab;    
    nrw = size(notiva,1);  
    newSS = [];
%     for iva = 1:nrw
    while iva <= nrw
        flag = 0;
        tiva = notiva(iva,:);
%         notiva(iva,:) = [];        
        tiva = repmat(tiva,size(notiva,1),1);
        if ~isempty(tiva)
        if notiva(iva,1) > 1 && notiva(iva,1) < 10
            tfr = abs(notiva(:,1)-tiva(:,1))<1e0;
            if any(tfr) && length(find(tfr)) > 1                  
                newSS = [newSS;tiva(1,1) ...
                             sum(notiva(abs(notiva(:,1)-tiva(:,1))<1e-2,2)) ...
                             sum(notiva(abs(notiva(:,1)-tiva(:,1))<1e-2,2))*100/ncnt];
                notiva(abs(notiva(:,1)-tiva(:,1))<1e-2,:) = [];                    
                iva = 1;
                flag = 1;
            elseif find(tfr) == 1
                newSS = [newSS;tiva(1,:)];
                notiva(tfr,:) = [];
                iva = 1;
                flag = 1;
            end 
        elseif notiva(iva,1) > 10
            tfr = abs(notiva(:,1)-tiva(:,1))<1e1;
            if any(tfr) && length(find(tfr)) > 1                  
                newSS = [newSS;tiva(1,1) ...
                               sum(notiva(tfr,2)) ...
                               sum(notiva(tfr,2))*100/ncnt];
                notiva(tfr,:) = [];                    
                iva = 1;
                flag = 1;
            elseif find(tfr) == 1
                newSS = [newSS;tiva(1,:)];
                notiva(tfr,:) = [];
                iva = 1;
                flag = 1;
            end 
        elseif notiva(iva,:) > 1e-3
            tfr = abs(notiva(:,1)-tiva(:,1))<1e-4;
            if any(tfr) && length(find(tfr)) > 1                  
                newSS = [newSS;tiva(1,1) ...
                               sum(notiva(tfr,2)) ...
                               sum(notiva(tfr,2))*100/ncnt];
                notiva(tfr,:) = [];                    
                iva = 1;
                flag = 1;
            elseif find(tfr) == 1
                newSS = [newSS;tiva(1,:)];
                notiva(tfr,:) = [];
                iva = 1;
                flag = 1;
            end 
        else
            tfr = abs(notiva(:,1)-tiva(:,1))<1e-6;
            if any(tfr) && length(find(tfr)) > 1                  
                newSS = [newSS;tiva(1,1) ...
                             sum(notiva(abs(notiva(:,1)-tiva(:,1))<1e-6,2)) ...
                             sum(notiva(abs(notiva(:,1)-tiva(:,1))<1e-6,2))*100/ncnt];
                notiva(abs(notiva(:,1)-tiva(:,1))<1e-6,:) = [];                    
                iva = 1;
                flag = 1;
            elseif find(tfr) == 1
                newSS = [newSS;tiva(1,:)];
                notiva(tfr,:) = [];
                iva = 1;
                flag = 1;
            end 
        end
        end
        nrw = size(notiva,1);
        if ~flag
            iva = iva + 1;
        end
    end
    multiSS{iv} = newSS;
%     display('iva complete');    
end




function [srate] = splitReg(trnmodel,ngene)

% function [Smodel,REffect,RConc,RLogic,srate] = splitReg(trnmodel,ngene)
%October 2013

%bindaffcell = cell(ngene,1);
srate = cell(ngene,1);
for igene = 1:ngene
    %bindaffcell{igene} = zeros(length(find(trnmodel.RS(igene,:))),1);
    srate{igene} = zeros(length(find(trnmodel.RS(igene,:))),1);
    tfcount = 1;
    nregterm = length(trnmodel.GeneRules{igene});
    iregterm = 1;
    newterms = {};    
    while iregterm <= nregterm
        [terms,xterms] = regexp(trnmodel.GeneRules{igene}{iregterm},'(\w+.?)(\W+?)+','tokens','split');
        if ~isemptyr(xterms)
            terms = [terms,cell(1,1)];
            terms{end}{1} = xterms{end};
        end
        newterms = [newterms,terms];
        iregterm = iregterm + 1;
    end
    %terms = newterms;
    nterms = length(newterms);
    protemp = cell(nterms,1);
    lterms = 1;
    
    while lterms <= nterms
        protemp{lterms} = newterms{lterms}{1};
        lterms = lterms + 1;
    end
    
    if nterms ~= length(find(trnmodel.Coefficient(igene,:)))
        nonprot = setdiff(trnmodel.Protein(logical(trnmodel.Coefficient(igene,:))),protemp);
        jterm = 1;
        while jterm <= length(nonprot)
            newterms{end}{2} = '|';
            newterms = [newterms,cell(1,1)];
            protindx = strcmp(nonprot{jterm},trnmodel.Protein);
            newterms{end}{1} = trnmodel.Protein{protindx};                     
            jterm = jterm + 1;
        end
    end
    nterms = length(newterms);
    kterm = 1;
    while kterm <= nterms
        tf = strcmp(newterms{kterm}{1},trnmodel.Protein);
        if any(tf)
            %bindaffcell{igene}(tfcount) = trnmodel.Coefficient(igene,tf);
            %bindaffcell{igene}(tfcount) = trnmodel.ProtConc(tf)/trnmodel.Coefficient(igene,tf);
            srate{igene}(tfcount) = trnmodel.srate(igene,tf);
            tfcount = tfcount + 1;
        end
        kterm = kterm + 1;
    end  
end
%====================== Original splitReg.m ===============================
% if ngene > 1
%     Smodel = struct();
%     Smodel.REffect = cell(ngene,1);
%     Smodel.RConc = cell(ngene,1);
%     Smodel.RLogic = cell(ngene,1);
%     %Smodel.RCoeff = cell(ngene,1);
%     Smodel.Srate = cell(ngene,1);
%     REffect = {};
%     RConc = {};
%     RLogic = {};
%     %RCoeff = {};    
%     %Srate = {};
%     flag = 1;
% else
%     Smodel = struct([]);
%     REffect = cell(ngene,1);
%     RConc = cell(ngene,1);
%     RLogic = cell(ngene,1);
%     %RCoeff = cell(ngene,1);
%     %Srate = cell(ngene,1);
%     flag = 0;
% end

% for igene = 1:ngene
%     %effect = trnmodel.RS(igene,logical(trnmodel.RS(igene,:)));
%     nreg = length(find(trnmodel.RS(igene,:)));
%     ireg = 1;
%     rconc = zeros(nreg,1);%To be replaced with intial regulatory protein concentrations
%     rlogic = zeros(nreg,1);
%     %rcoeff = zeros(nreg,1);
%     reffect = zeros(nreg,1);
%     %srate = zeros(nreg,1);
%     while ireg <= nreg
%         reffect(ireg) = trnmodel.RS(igene,trnmodel.Order(igene,:)==ireg);
%         rconc(ireg) = trnmodel.ProtConc(trnmodel.Order(igene,:)==ireg);
%         rlogic(ireg) = trnmodel.Logic(igene,trnmodel.Order(igene,:)==ireg);
%         %rcoeff(ireg) = trnmodel.Coefficient(igene,trnmodel.Order(igene,:)==ireg);
%         %srate(ireg) = trnmodel.srate(igene,trnmodel.Order(igene,:)==ireg);               
%         ireg = ireg + 1;
%     end
%     if flag
%         Smodel.REffect{igene} = reffect;
%         Smodel.RConc{igene} = rconc;
%         Smodel.RLogic{igene} = rlogic;
%         %Smodel.RCoeff{igene} = rcoeff;
%         %Smodel.Srate{igene} = srate;
%     else
%         REffect{igene} = reffect;
%         RConc{igene} = rconc;
%         RLogic{igene} = rlogic;
%         %RCoeff{igene} = rcoeff;
%         %Srate{igene} = srate;
%         
%     end        
%     %dbstop in pactivity.m
%     %promoteract(igene) = pactivity(nreg,rconc,rcoeff,reffect,rlogic);
% end
end
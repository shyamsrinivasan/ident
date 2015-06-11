function [bindaffcell] = getbindingaffinity(trnmodel,ngene,protconc) 
%Add it to splitReg.m??
%The end result of both these functions seems to be the same
%However calculation of binding affinity is to be done at every iteration
%of the ODE numerical solver. Hence better to write it as a separate
%function
%Getting srate everytime getbindingaffinity is called is an additional
%computational burden. Consider an alternative like splitReg.m/ but not
%associated with order
%October 2013
bindaffcell = cell(ngene,1);
%srate = cell(ngene,1);
for igene = 1:ngene
    bindaffcell{igene} = zeros(length(find(trnmodel.RS(igene,:))),1);
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
    %nterms = length(newterms);
    protemp = cell(length(newterms),1);
    lterms = 1;
    
    while lterms <= length(newterms)
        protemp{lterms} = newterms{lterms}{1};
        lterms = lterms + 1;
    end
    
    if length(newterms) ~= length(find(trnmodel.Coefficient(igene,:)))
        nonprot = setdiff(trnmodel.Protein(logical(trnmodel.Coefficient(igene,:))),protemp);
        jterm = length(nonprot);
        while jterm >= 1
            newterms{end}{2} = '|';
            newterms = [newterms,cell(1,1)];
            protindx = strcmp(nonprot{jterm},trnmodel.Protein);
            newterms{end}{1} = trnmodel.Protein{protindx};                     
            jterm = jterm - 1;
        end
    end
    nterms = length(newterms);
    kterm = 1;
    while kterm <= nterms
        tf = strcmp(newterms{kterm}{1},trnmodel.Protein);
        if any(tf)
            %bindaffcell{igene}(tfcount) = trnmodel.Coefficient(igene,tf);
            %bindaffcell{igene}(tfcount) = protconc(tf)/trnmodel.Coefficient(igene,tf);
            bindaffcell{igene}(kterm) = protconc(tf)/trnmodel.Coefficient(igene,tf);
            %srate{igene}(tfcount) = trnmodel.srate(igene,tf);
            %tfcount = tfcount + 1;
        end
        kterm = kterm + 1;
    end  
end
end

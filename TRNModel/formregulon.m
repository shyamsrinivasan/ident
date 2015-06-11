function [regulon] = formregulon(trnmodel,igene,nregs)
regulon = cell(1,1);
regs = trnmodel.Protein(logical(trnmodel.RS(igene,:)));
rmat = trnmodel.RS(igene,logical(trnmodel.RS(igene,:)));
ireg = 1;
while ireg <= nregs
    switch rmat(ireg)
        case 1
            regulon{1} = [regulon{1},sprintf('%s(+)',regs{ireg})];
        case -1
            regulon{1} = [regulon{1},sprintf('%s(-)',regs{ireg})];
        case 2
            regulon{1} = [regulon{1},sprintf('%s(+/-)',regs{ireg})];
    end
    ireg = ireg + 1;
    %newC{igene} = regulon{1};
end 
end



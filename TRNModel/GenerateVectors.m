function [TF,DNA,Gene1,PGenes,EGene,EPSGene]=GenerateVectors(model,pgenes)
%Removed non-essential commented statements April 25th

PGenes = pgenes;

%Model has the initial E matrix
%ETF = model.ETF;
%EDNA1 = model.EDNA;
EGene = model.EGene;
EPSGene = model.EPSGene;

nofg = size(EGene,1);
nofinputs = size(PGenes,1);


TF = ones(nofg,1);
Gene1 = ones(size(model.GeneName,1),1);
DNA = ones(nofg,1);

end
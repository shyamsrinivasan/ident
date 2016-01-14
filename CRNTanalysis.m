function [SBMLmodel,result] = CRNTanalysis(FBAmodel)

SBMLmodel = converttoSBMLformat(FBAmodel);

%CRNT analysis
result = model_analysis(SBMLmodel);

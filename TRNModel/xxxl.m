%Extracting Activation & Repression Coefficients in One go
function [trnmodel,ngene,tnreg,tnact,tnrep,tnactrep] = Tmodel(rxfname,regtype)

fileid = fopen(rxfname);

if fileid == -1
    fprintf('File %s cannot be opened.', rxfname);
    trnmodel = struct([]);
    return;
end

C = textscan(fileid, '%s%s%s%f%f', 'Delimiter', '\t', 'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);

trnmodel.Gene = C{1};
ngene = length(trnmodel.Gene);

sym = {'(+)';'(-)';'(+/-)'}; %To be tied-up with regtype in the function input arguments

trnmodel.Activator = cell(ngene,1);
trnmodel.Repressor = cell(ngene,1);
trnmodel.ActRep = cell(ngene,1);

trnmodel.ActiveCoeff = cell(ngene,1);
trnmodel.RepressCoeff = cell(ngene,1);
trnmodel.ActRepCoeff = cell(ngene,1);

%AcCoeff = cell(ngene,1);
%RepCoeff = cell(ngene,1);
%OthCoeff = cell(ngene,1);
tnreg = 0;
tnact = 0;
tnrep = 0;
tnactrep = 0;
%tncoeff = 0;
pflag = 0;
for igene = 1:ngene
    [terms,startindx,endindx] = regexp(C{2}{igene},'(\w+.?)(\(\W+.?\))+','tokens');
    [coeff,startindx2,endindx2] = regexp(C{3}{igene},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
    
    nterms = length(terms);
    ncoeff = length(coeff);
    %flag = 0; Consistency Check Flags - To be used later
        
    if ncoeff ~= nterms 
        fprintf('Number of Coefficients ~= Number of Regulators \n');
        fprintf('Gene %s will have default values \n',C{1}{igene});
        if ~pflag
            fprintf('Switching to Default Values \n');
            fprintf('Default Values \n Activation Coeff. %d\n Repression Coeff. %d\n Other Coeff. %d\n',...
                defaccoeff,defrepcoeff,defacrepcoeff);
            pflag = 1;
        end
        
        Regulators = C{2}{igene};
        if ~iscell(Regulators)
            Regulators = cellstr(Regulators);
        end
        Coeff = initRegCoeff(Regulators,length(Regulators));
        [coeff,startindx2,endindx2] = regexp(Coeff{1},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
        %flag = 1; Consistency Check Flags - To be used later
    end
    
    iterm = 1;
    actindx = 0;
    repindx = 0;
    acrepindx = 0;
    activator = cell(nterms,1);
    activecoeff = zeros(nterms,1);
    repressor = cell(nterms,1);
    represscoeff = zeros(nterms,1);
    actrep = cell(nterms,1);
    actrepcoeff = zeros(nterms,1);
     
    while iterm <= nterms %%&& flag == 0
        %regulon = terms{iterm}{1};
        %effect = terms{iterm}{2};
        if any(strcmp(terms{iterm}{2},sym{1})) %Activation (+)
            actindx = actindx + 1;
            activator{actindx} = terms{iterm}{1};
            activecoeff(actindx) = str2double(coeff{iterm}{1});
            tnreg = tnreg + 1;
            %tncoeff = tncoeff + 1;
            tnact = tnact + 1;
        elseif any(strcmp(terms{iterm}{2},sym{2})) %Repression (-)
            repindx = repindx + 1;
            repressor{repindx} = terms{iterm}{1};
            represscoeff(repindx) = str2double(coeff{iterm}{1});
            tnreg = tnreg + 1;
            %tncoeff = tncoeff + 1;
            tnrep = tnrep + 1;
        elseif any(strcmp(terms{iterm}{2},sym{3})) %Activation/Repression (+/-)
            acrepindx = acrepindx + 1;
            actrep{acrepindx} = terms{iterm}{1};
            actrepcoeff(acrepindx) = str2double(coeff{iterm}{1});
            tnreg = tnreg + 1;
            %tncoeff = tncoeff + 1;
            tnactrep = tnactrep + 1;
        end
        iterm = iterm + 1;
    end
      
    activator = activator(~cellfun('isempty',activator));
    repressor = repressor(~cellfun('isempty',repressor));
    actrep = actrep(~cellfun('isempty',actrep));
    
    trnmodel.Activator{igene} = activator;
    trnmodel.Repressor{igene} = repressor;
    trnmodel.ActRep{igene} = actrep;
    
    activecoeff(activecoeff == 0) = [];
    represscoeff(represscoeff == 0) = [];
    actrepcoeff(actrepcoeff == 0) = [];
    
    trnmodel.ActiveCoeff{igene} = activecoeff;
    trnmodel.RepressCoeff{igene} = represscoeff;
    trnmodel.ActRepCoeff{igene} = actrepcoeff;
    
end

end
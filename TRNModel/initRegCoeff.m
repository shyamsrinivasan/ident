%function Coeff = initRegCoeff(Regulon,ngene,defparval,dbfname,sheetnum)

%Writing Activation & Repression Coefficients for unknowns/undefined
%******Inputs
%Regulon -   String of regulators for genes or cell array of genes
%ngene -     Number of genes in Regulon
%******Optional Inputs
%defparval - Structure of default parameter values for TRN models
%dbfname -   If values are to be written to an Excel file, specifies 
%            the file name
%sheetnum -  Sheet number on the Excel file specified using dbfname
%******Outputs
%Coeff -     Cell array of coefficients for regulators in Regulon and the 
%            symbol for their corresponding effect as specified through 
%            Regulon
%October 2013 - version 1.0 
%July 2014 - Superseded by coeffinit v2.0  
%**************************************************************************             
function Coeff = initRegCoeff(Regulon,ngene,defparval,dbfname,sheetnum)

%Still to correct the write to file portion
%December 05 2013

if nargin < 5
    sheetnum = 'Sheet1';
end
if nargin < 4
    %fprintf('No filename given \n Not writing results to excel file \n');
    wflag = 0;
end

if nargin < 3
    defparval.accoeff = 1e-9;
    defparval.repcoeff = 1e-9;
    defparval.rephill = 2;
    defparval.dualcoeff = 1e-9;
    defparval.brate = 1.66e-6;
    defparval.srate = 1.66e-6;
    defparval.trate = 4.94e-6;
    defparval.drate = 1.44e-3;
    defparval.kmax = 5e-6;
    defparval.ks = 1e-6;
    defparval.pdrate = 1.44e-3;
end
if nargin < 2
    ngene = length(Regulon);
end
    
sym = {'(+)';'(-)';'(+/-)'};
Coeff = cell(ngene,1);

for igene = 1:ngene
    if ~isempty(Regulon{igene})
        terms = regexp(Regulon{igene},'(\w+.?)(\(\W+\))+','tokens');
        nterms = length(terms);
        Coeff = cell(1,nterms);
        iterm = 1;
        while iterm <= nterms 
            %regulon = terms{iterm}{1};
            %effect = terms{iterm}{2};
            if any(strcmp(terms{iterm}{2},sym{1})) %Activation
                Coeff{iterm}{1} = defparval.accoeff;
                Coeff{iterm}{2} = '(+)';
                %Coeff{igene} = [Coeff{igene},sprintf('%g(+)',defaccoeff)];
                %tnact = tnact + 1;
            elseif any(strcmp(terms{iterm}{2},sym{2})) %Repression
                Coeff{iterm}{1} = defparval.repcoeff;
                Coeff{iterm}{2} = '(-)';
                %Coeff{igene} = [Coeff{igene},sprintf('%g(-)',defrepcoeff)];
                %tnrep = tnrep + 1;
            elseif any(strcmp(terms{iterm}{2},sym{3})) %Activation/Repression
                Coeff{iterm}{1} = defparval.accoeff;
                Coeff{iterm}{2} = '(+/-)';
                %Coeff{igene} = [Coeff{igene},sprintf('%g(+/-)',defacrepcoeff)];
                %tnactrep = tnactrep + 1;
            end
            iterm = iterm + 1;
        end
        if any(wflag)
%             cellindx = sprintf('C%d',igene+1);
%             xlswrite(dbfname,cellstr(Coeff{igene}),sheetnum,cellindx);
        end
    end
end
end
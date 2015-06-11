% function [strates,btrate,pflag] = srateinit(trnmodel,Regulators,defparval,pflag,wflag,dbfname,sheetnum)
%**************************************************************************
%Function to initialize stimulated transcription rates for regulators for a
%single gene
%OCtober 2013 - version 1.0

%**************************************************************************
function [strates,btrate,pflag] = srateinit(trnmodel,Regulators,defparval,pflag,wflag,dbfname,sheetnum)
if nargin < 7
    sheetnum = 'Sheet1';
end
if nargin < 6
    %fprintf('No filename given \n Not writing results to excel file \n');
    wflag = 0;
end
if nargin < 5
    wflag = 0;
end
if nargin < 4
    pflag = 0;
end
if nargin < 3
    defstrate = 0.000023;%Was previously 0.16 
    defbtrate = 2.3e-6;%Was previously 0.0133 
else
    defstrate = defparval.srate;
    defbtrate = defparval.brate;
end


if ~pflag
    fprintf('Default Rates\n Transcription:%4.5g\n Translation:%4.5g\n Decay:%4.5g\n',...
        defparval.srate,defparval.trate,defparval.drate);
    pflag = 1;
end
    
%Function is supposed to be called within a loop and the for loop within
%the function should be done away with as it cannot function as intended
%Nov 25 2013

ngene = length(Regulators);
%strates = cell(ngene,1);
%btrate = zeros(ngene,1);
%srterms = cell(ngene,1);
%transrate = cell(ngene,1);

%This code only supports call for ngene = 1!
for igene = 1:ngene %ngene > 1
    %transrate{igene} = cell(2,1);
    btrate = defbtrate;
%     if ~isempty(Regulators) %ngene = 1
    if ~isempty(Regulators{igene}) %ngene > 1
        %%ngene > 1
        terms = regexp(Regulators{igene},'(\w+.?)(\(\W+\))+','tokens');
        %ngene = 1
%         terms = regexp(Regulators{1},'(\w+.?)(\(\W+\))+','tokens');

        nterms = length(terms);
        
        %Possible replacement for the statements below to determine srterms
        %Nov 25 2013
        
        strates = cell(1,nterms);
        if wflag
%             srterms = [];
%             srterms = char(srterms);
            regterms = [];
            regterms = char(regterms);
        end
               
        iterm = 1;
        while iterm <= nterms 
            switch terms{iterm}{2}%Possibility of introducing different 
                                  %rates for activation and inhibition
                case '(+)'
                    %btrates{igene} = [btrates{igene},sprintf('%g(+)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%8.6f(+)',defstrate)];
                    strates{iterm} = {mat2str(defstrate);'(+)'};
                    if wflag
                        
                        %srterms = [srterms,sprintf('%5.4g(+)',defparval.accoeff)];
%                         srterms = [srterms,sprintf('%5.4g(+)',defparval.srate)];
%                         regterms = [regterms,sprintf('%s(+)',terms{iterm}{1})];
                    end
                    
                case '(-)'
                    %btrates{igene} = [btrates{igene},sprintf('%g(-)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%8.6f(-)',defstrate)];
                    strates{iterm} = {mat2str(defstrate);'(-)'};
                    if wflag
                        %srterms = [srterms,sprintf('%5.4g(-)',defparval.repcoeff)];
                        %srterms = [srterms,sprintf('%5.4g(-)',defparval.srate)];
%                         regterms = [regterms,sprintf('%s(-)',terms{iterm}{1})];
                    end
                    
                case '(+/-)'
                    %btrates{igene} = [btrates{igene},sprintf('%g(+/-)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%8.6f(+/-)',defstrate)];
                    strates{iterm} = {mat2str(defstrate);'(+/-)'};
                    if wflag
                        %srterms = [srterms,sprintf('%5.4g(+/-)',defparval.dualcoeff)];
%                         srterms = [srterms,sprintf('%5.4g(+/-)',defparval.srate)];
%                         regterms = [regterms,sprintf('%s(+/-)',terms{iterm}{1})];
                    end                    
            end
            iterm = iterm + 1;
        end
        
        %File writing operating needs correction to incorporate above changes
        %Removed until further notice!
        
        if wflag 
%             cellindx1 = sprintf('C%d',igene+1);
%             xlswrite(dbfname,cellstr(btrates{igene}),sheetnum,cellindx1);
            fprintf('Writing Data for Gene %d\n',igene);
            cellindx2 = sprintf('A%d',igene);
            xlswrite(dbfname,cellstr(srterms),sheetnum,cellindx2);
        end
    end
end
    %[brterms] = regexp(btrates{igene},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
    %[srterms] = regexp(strates{igene},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
    
    %Remove above statement and replace with the following statement
    %Nov 25 2013
    
    %srterms = strates{igene};
    %srterms = strates;
    
    %srterms{igene} = srterms1;
    
%end  


   
end
function [brterms,pflag] = brateinit(trnmodel,Regulators,defbtrate,defstrate,wflag,dbfname,sheetnum,pflag)
if nargin < 8
    pflag = 1;
end
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
% if nargin < 4
%     defstrate = 0.16;
% end
if nargin < 3
    defbtrate = 0.0133;
end

ngene = length(Regulators);
%strates = cell(ngene,1);
btrates = cell(ngene,1);
%transrate = cell(ngene,1);
for igene = 1:ngene
    %transrate{igene} = cell(2,1);
    if ~isempty(Regulators{igene})
        terms = regexp(Regulators{igene},'(\w+.?)(\(\W+\))+','tokens');
        nterms = length(terms);
        
        iterm = 1;
        while iterm <= nterms 
            switch terms{iterm}{2}
                case '(+)'
                    btrates{igene} = [btrates{igene},sprintf('%g(+)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%g(+)',defstrate)];
                case '(-)'
                    btrates{igene} = [btrates{igene},sprintf('%g(-)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%g(-)',defstrate)];
                case '(+/-)'
                    btrates{igene} = [btrates{igene},sprintf('%g(+/-)',defbtrate)];
                    %strates{igene} = [strates{igene},sprintf('%g(+/-)',defstrate)];
            end
            iterm = iterm + 1;
        end
        if wflag
            cellindx1 = sprintf('C%d',igene+1);
            xlswrite(dbfname,cellstr(btrates{igene}),sheetnum,cellindx1);
%             cellindx2 = sprintf('D%d',igene+1);
%             xlswrite(dbfname,cellstr(strates{igene}),sheetnum,cellindx2);
        end
    end
    [brterms] = regexp(btrates{igene},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
    %[srterms] = regexp(strates{igene},'(\d*\.?\d+)(\(\W+.?\))+','tokens');
    
%     transrate{igene}{1} = brterms;
%     transrate{igene}{2} = srterms;
end
    
    
    
    
 
    


   
end
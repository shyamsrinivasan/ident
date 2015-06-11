%Get only Final Steady state values from the dynamic solution set.
function [finalSSdata,finalTdata,f_names] = getssval(Solution,ngene,tnreg)
% Yet to Finalized Below!!
% 
f_names = {'initSS'};
%if exist('Solution','var')
names = fieldnames(Solution);
nfields = length(names);
if nfields > 1
    for if_name =  1:nfields-1
        %newf_name = sprintf('pertb%d_SS',if_name);
        f_names = [f_names;sprintf('pertb%d',if_name)];
    end
end

finalSSdata = zeros(ngene+tnreg,nfields);
finalTdata = zeros(ngene+tnreg,nfields);
%finalSSdata1 = [];

SStime = zeros(ngene+tnreg,nfields);
for if_name = 1:nfields
    if any(strcmp(f_names{if_name},names))
        fieldname = names{strcmp(f_names{if_name},names)};
        finalSSdata(:,if_name) = Solution.(fieldname){2}(:,end);
        finalTdata(:,if_name) = Solution.(fieldname){1}(end);
        
        %finalSSdata1 = [finalSSdata1,Solution.(fieldname){2}(:,15211)];
        
        %Determining time at which SS is achieved
        
        %             ntimepts = length(Solution.(fieldname){1});
        %             errtol = zeros(ngene+tnreg,1);
        %             errtol = 10^floor(log10(Solution.(fieldname){2}(:,1)));
        
        %             for itpts = 2:ntimepts
        %                 for ivar = 1:ngene+tnreg
        %                     if abs(Solution.(fieldname){2}(ivar,itpts-1) -...
        %                             Solution.(fieldname){2}(ivar,itpts)) <= 1e-6
        %                         SStime(ivar,if_name) = Solution.(fieldname){1}(ivar,itpts);
        %                     end
        %                 end
        %             end
        
        %                     SStime(:,if_name) = Solution.(fieldname){1}(
        %
        %
        
    end
end
%end
end
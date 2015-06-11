%function [AllDynSolution] = storeSolution(Solution,AllDynSolution,...
%                            trnmodel,ngene,gname,pname)
%**************************************************************************
%Store all dynamic information from multiple simulation experiments
%*****Inputs
%Solution           Solution structure (in the form of output given by
%                   trnperturbation.m
%AllDynSolution     Structure of previously stored dynamic information during
%                   a previous call to storeSolution.Empty struture with no
%                   fields otherwise
%trnmodel
%ngene
%gname              Cell array of gene names whose dynamic profiles are to 
%                   be stored
%pname              Cell array of protein names whose dynamic profiles are 
%                   to be stored
%*****Outputs
%AllDynSolution     Structure of dyanmic information on a gene by gene 
%                   basis for all genes and proteins in gname and pname 
%                   respectively
%**************************************************************************
%Notes:
%January 9th 2014
%Need to add support for multiple gene names (> 1)
%Need to add support for protein profiles
%January 11th 2014
%Support added for multiple genes as well as protein names to be done in a
%single call to storeSolution
%Should be called within a loop to store results from multiple simulations
%Limitation: Can only store data on a gene by gene basis
function [AllDynSolution,AllSSsolution] = storeSolution(Solution,...
                AllSSsolution,AllDynSolution,trnmodel,ngene,gname,pname)

if nargin < 2
    AllDynSolution = struct();
end

if ~isempty(Solution)
    sfnames = fieldnames(Solution);
    nfields = length(sfnames);
    for jfield = 1:nfields
        tss = Solution.(sfnames{jfield}){1};
        gss = Solution.(sfnames{jfield}){2}(1:ngene,:)';
        pss = Solution.(sfnames{jfield}){2}(ngene+1:end,:)';
        
        %Update final SS solution
        %Store the final SS values in a separate structure AllSSsolution
        
        
        if isfield(AllSSsolution,'finalT')
            if isfield(AllSSsolution.finalT,sfnames{jfield})
                time_var = AllSSsolution.finalT.(sfnames{jfield});
                time_var = [time_var;tss(end)];
                AllSSsolution.finalT.(sfnames{jfield}) = [];
                AllSSsolution.finalT.(sfnames{jfield}) = time_var;
            else
                AllSSsolution.finalT.(sfnames{jfield}) = tss(end);
            end
        else
            AllSSsolution.finalT.(sfnames{jfield}) = tss(end);
        end
        
        if isfield(AllSSsolution,'finalSS')
            if isfield(AllSSsolution.finalSS,sfnames{jfield})
                exp_var = AllSSsolution.finalSS.(sfnames{jfield});
                AllSSsolution.finalSS.(sfnames{jfield}) = [];
                exp_var = [exp_var,Solution.(sfnames{jfield}){2}(:,end)];
                AllSSsolution.finalSS.(sfnames{jfield}) = exp_var; 
            else
                AllSSsolution.finalSS.(sfnames{jfield}) =...
                    Solution.(sfnames{jfield}){2}(:,end);
            end
        else
            AllSSsolution.finalSS.(sfnames{jfield}) =...
                Solution.(sfnames{jfield}){2}(:,end);
        end
            
                
%             ssfname = fieldnames(AllSSsolution.finalT);
%             time_var = AllSSsolution.finalT;
%             AllSSsolution.finalT = [];
%             time_var = [time_var;tss(end)];
%             AllSSsolution.finalT = time_var;
%         else
%             AllSSsolution.finalT = tss(end);
%         end
%         if isfield(       
            
            
        
%         updateSSsolution(tss,Solution.(sfnames{jfield}{2}));
        
        if ~isempty(gname)
            for igene = 1:length(gname)
                tfg = strcmp(gname{igene},trnmodel.Gene);
                if any(tfg)
                updateSolutionStrutcure(trnmodel.Gene{tfg},...
                                        sfnames{jfield},tss,gss,tfg);
                end
            end
        end
        if ~isempty(pname)
            for iprot = 1:length(pname)
                tfp = strcmp(pname{iprot},trnmodel.Protein);
                if any(tfp)
                updateSolutionStrutcure(trnmodel.Protein{tfp},...
                                        sfnames{jfield},tss,pss,tfp);
                end
            end
        end
    end
else
    fprintf('No Solution to Save. Solution Structure is empty\n');
end
% function updateSSsolution(time,expression)
%     temp_var = finalSSsolution;
%     temp_var = [temp_var, expression(:,end)];
%     time_var = finalT;
%     finalT = [];
%     time_var = [time_var;time(end)];
%     finalT = time_var;
%     %Algorithm
%     %The final SS values are in the last column/rows only
%     %SStime = [SStime,time(end)];
%     %expression = [expression,expression(:,end)];
%     
% end

function updateSolutionStrutcure(gp_name,field_name,timevec,expvec,indx)
    adsfnames = fieldnames(AllDynSolution);
    adsfndx = strcmp(sprintf('%s_%s',gp_name,field_name),adsfnames);
    if any(adsfndx)%Field Exists
        temp_var = AllDynSolution.(adsfnames{adsfndx});
        AllDynSolution.(adsfnames{adsfndx}) = [];
        temp_var = [temp_var timevec expvec(:,indx)];
        AllDynSolution.(adsfnames{adsfndx}) = temp_var;
        
        %flag = 1;%Do not do anything or store values
    else%Create such a field & allocate memory for the future
        temp_var = [timevec expvec(:,indx)];
        AllDynSolution.(sprintf('%s_%s',gp_name,field_name)) = temp_var;
    end   
    
    
end



end
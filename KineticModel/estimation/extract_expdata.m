function [exp_xss,exp_fss] = extract_expdata(exp_sol,id,nval,nc,nf,npert)


% experimental data
if isempty(id)
    if ~isempty(nc)
        exp_xss = zeros(nval*nc,npert);
        if iscell(exp_sol)
            for icell = 1:nval
                if isfield(exp_sol{icell},'xss')
                    exp_xss((icell-1)*nc+1:nc*icell,:) = cat(2,exp_sol{icell}.xss);
                end
            end
        else
            if isfield(exp_sol,'xss')
                exp_xss = cat(2,exp_sol.xss);
            end    
        end
    end
    
    if ~isempty(nf)
        exp_fss = zeros(nval*nf,npert);
        if iscell(exp_sol)
            for icell = 1:nval
                if isfield(exp_sol{icell},'fss')
                    exp_fss((icell-1)*nf+1:nf*icell) = cat(2,exp_sol{icell}.fss);
                end
            end
        else
            if isfield(exp_sol,'fss')
                exp_fss = cat(2,exp_sol.fss);
            end
        end
    end
else
    if iscell(exp_sol)
        exp_xss = cat(2,exp_sol{id}.xss);
        exp_fss = cat(2,exp_sol{id}.fss);
    else
        exp_xss = cat(2,exp_sol.xss);
        exp_fss = cat(2,exp_sol.fss);
    end
end
%Extract Data either on a Expt by Expt basis or an entire set of
%Experiments
function [varargout] = ExtractData(call,varargin)
if nargin < 5
    varargin{4} = {};
end
if call == 1
    %Assign Metabolite concentrations
elseif call == 2
    %Assign varname for plotting => MetabConc = {}
    mName = varargin{1};
    if ~isempty(mName)
        varname = parsestring(mName);
    end
    varargout{1} = varname;
elseif call == 3
    %Assign parameters?
    if isstruct(varargin{1})
        model = varargin{1};
        %Extract parameters
        parname = parsestring(varargin{2});
        parval = parsestring(varargin{3});
        if ~isempty(varargin{4})
            parapp1 = parsestring(varargin{4});
        else
            parapp1 = {};
        end
        napp1 = length(parapp1);
        gene = cell(napp1,1);
        prot = cell(napp1,1);
        for iapp1 = 1:napp1
            hpos = strfind(parapp1{iapp1},'-');
            gene{iapp1} = parapp1{iapp1}(1:hpos-1);
            prot{iapp1} = parapp1{iapp1}(hpos+1:end);
        end
        for ipar = 1:length(parname)
            if strcmpi(parname{ipar},'ac')%activation 
                [model] = change_pmeter(parval{ipar},gene,prot,napp1,model);                
            elseif strcmpi(parname{ipar},'ic')%inhibition 
                [model] = change_pmeter(parval{ipar},gene,prot,napp1,model);
            elseif strcmpi(parname{ipar},'brate')%basal
                [model] = change_pmeter(parval{ipar},gene,prot,napp1,model);
            elseif strcmpi(parname{ipar},'hc')%hill
            else%kcat,Vuptake,Vefflux,gmax
                model.(parname{ipar}) = str2double(parval{ipar});
                %setting default values
                if ~isfield(model,'Vuptake')
                    [~,Vuptake] = find(model.S(:,model.Vic_exind)>0);
                    Vuptake = model.Vic_exind(Vuptake);
                    model.Vuptake = zeros(length(Vuptake),1);
                end
                if ~isfield(model,'Vefflux')
                    [~,Vexcrt] = find(model.S(:,Vic_exind)<0);
                    Vexcrt = Vic_exind(Vexcrt);
                    model.Vefflux = zeros(length(Vexcrt),1);
                end
                if ~isfield(model,'kcat')
                    model.kcat = 10000;
                end
                if ~isfield(model,'gmax')
                    model.gmax = 0.8;
                end
            end
        end
    end  
    varargout{1} = model; 
    varargout{2} = parname;
    varargout{3} = parval;
elseif call == 4
    %Sample initial concentrations
    varargout{1} = str2double(varargin{1});
end

function [model] = change_pmeter(spar_val,gene_vec,prot_vec,napp1,model)
[Coefficient,brate] = parameter_return(model.allpar,model);
for iapp1 = 1:napp1
    if ~isempty(gene_vec)
        tfg = strcmpi(gene_vec{iapp1},model.Gene);
    else
        tfg = 0;
    end
    if ~isempty(prot_vec)
        tfp = strcmpi(prot_vec{iapp1},model.Regulators);
    else
        tfp = 0;
    end
    if any(tfg) && any(tfp)%for ac and ic
        if model.RS(tfg,tfp) ~= 0
            if Coefficient(tfg,tfp) ~= 0
                Coefficient(tfg,tfp) = str2double(spar_val);
            end
        end
    elseif any(tfg)%for brate
        brate(tfg) = str2double(spar_val);
    end
end
newstruct.Coefficient = Coefficient;
newstruct.brate = brate;
model.allpar = parameter_vector(newstruct,length(model.Gene));
return

function [pstring] = parsestring(string)
compos = strfind(string,',');
nterm = length(compos)+1;
pstring = cell(nterm,1);
for iterm = 1:nterm
    if iterm == 1 && nterm > 1
        pstring{iterm} = string(1:compos(iterm)-1);
    elseif iterm > 1 && iterm < nterm 
        pstring{iterm} = string(compos(iterm-1)+1:compos(iterm)-1);
    elseif iterm == nterm && nterm > 1
        pstring{iterm} = string(compos(iterm-1)+1:end);
    else%iterm == nterm && nterm = 1
        pstring{iterm} = string(1:end);
    end
end
return
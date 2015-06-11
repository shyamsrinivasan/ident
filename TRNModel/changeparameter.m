% function [trnmodel] =
% changeparameter(trnmodel,metabprot,parstruct,genename,protname,regname,
%                 parname,newparval)
%**************************************************************************
%Change Parameters for specific/all genes or proteins in model
%*****Inputs
%GeneStruct - Structure with gene related para,eter fields  
%Fields
%gname      Gene name
%brate      Array of Basal mRNA transcription rate with size such that
%           number of rows of brate = length(GeneStruct.gname)
%drate      mRNA decay rate with size such that
%           number of rows of drate = length(GeneStruct.gname)
%allreg     Field to specify wether individual parameters are to be changed
%           or all parameters in the model are to be changed for regulators
%           corresponding to a gene
%allgene    Field to specify wether individual parameters are to be changed
%           or all parameters in the model are to be changed for Genes in
%           the model
%regulators Regulators associated with the gene
%srate      Stimulated transcription rate
%accoeff    Activation Coefficient
%repcoeff   Repression Coefficient
%dualcoeff  Dual regulation coefficient

%ProtStruct - Structure with protein related paarmeter field
%Fields
%pname      Protein Name
%trate      Protein tramslation rate
%pdrate     Protein decay rate
%allprot    Field to specify wether individual parameters are to be changed
%           or all parameters in the model are to be changed for proteins
%           in the models
%geneass    Genes associated with this protein in terms of regulation
%*****Outputs
%**************************************************************************
%Notes:
%December 13th 2013
%Need to support the ability to change one or all parameters in a given
%model
%January 11 2014
%Can change parameters srate and coeff for all genes if GeneStruct.gname is
%not available
%Still does not support change of protein parameters

% function [trnmodel] = changeparameter(trnmodel,metabprot,parstruct,
%                       genename,protname,regname,parname,newparval)
function [trnmodel] =...
    changeparameter(trnmodel,defparval,regprotein,GeneStruct,ProtStruct)
if nargin < 5
    ProtStruct = struct([]);
end
GeneFields = {'gname';'brate';'srate';'drate';'acccoeff';'repcoeff';'dualcoeff';...
              'regulators'};
ProtFields = {'pname';'trate';'pdrate'};  

if ~isempty(GeneStruct)
    if isfield(GeneStruct,'gname') && ~isfield(GeneStruct,'allgene')
        GeneStruct.allgene = 0;
        if isfield(GeneStruct,'regulators') && ~isfield(GeneStruct,'allreg')
            GeneStruct.allreg = zeros(length(GeneStruct.gname),1);
        else
            GeneStruct.allreg = ones(length(GeneStruct.gname),1);
        end            
    else
        GeneStruct.allgene = 1;
        if isfield(GeneStruct,'regulators')
            GeneStruct.allreg = 0;
        else
            GeneStruct.allreg = 1;
        end
        fprintf('No gene name given.');
        fprintf('Specified parameters will be changed for all genes\n');
    end
    
    if ~GeneStruct.allgene
        ngene = length(GeneStruct.gname);            
    else%GeneStruct.allgene = 1
        ngene = length(trnmodel.Gene);
        GeneStruct.gname = trnmodel.Gene;
        if GeneStruct.allgene
            GeneStruct.allreg = ones(length(GeneStruct.gname),1);
        else
            GeneStruct.allreg = zeros(length(GeneStruct.gname),1);
        end
            
    end
    
    for igene = 1:ngene
        tfg = strcmp(GeneStruct.gname{igene},trnmodel.Gene);
        if isfield(GeneStruct,'brate')
            if length(GeneStruct.gname) == length(GeneStruct.brate)
                ratechange({'brate'},tfg,GeneStruct.brate(igene));
            else
                %Mismatch between # of rates & # of genes. So use default
                %values for all genes
            end
        end
        if isfield(GeneStruct,'drate')
            if length(GeneStruct.gname) == length(GeneStruct.drate)
                ratechange({'drate'},tfg,GeneStruct.drate(igene));
            else
                %Mismatch between # of rates & # of genes. So use default
                %values for all genes
            end
        end
        if isfield(GeneStruct,'srate')
            if ~GeneStruct.allgene            
                ratechange({'srate'},tfg,...
                            GeneStruct.srate{igene},...
                            GeneStruct.regulators{igene},...
                            GeneStruct.allreg(igene));
            else
                ratechange({'srate'},tfg,...
                            [],...
                            trnmodel.Regulator{igene},...
                            GeneStruct.allreg(igene));
            end
        end
        if isfield(GeneStruct,'coeff')
            if ~GeneStruct.allgene            
                ratechange({'coeff'},tfg,...
                            GeneStruct.coeff{igene},...
                            GeneStruct.regulators{igene},...
                            GeneStruct.allreg(igene));
            else
                ratechange({'coeff'},tfg,...
                            [],...
                            trnmodel.Regulator{igene},...
                            GeneStruct.allreg(igene));
            end       
        end
        %             if isfield(GeneStruct,'repcoeff')
        %                 coeffchange(GeneFields(gfindx),GeneStruct{igene});
        %             end
        %             if isfield(GeneStruct,'dualcoeff')
        %                 coeffchange(GeneFields(gfindx),GeneStruct{igene});
        %             end
        
    end
    
end  

function flag = ratechange(par_name,indx,par_data,regulators,allreg)
    if nargin < 4
        regulators = {};
    end    
    
    par_valflag = 0;
    flag = 0;
       
    if strcmp(par_name,'coeff') || strcmp(par_name,'srate') 
        nreg = length(regulators);
        ireg = 1;
        while ireg <= nreg
            tfreg = strcmp(regulators{ireg},trnmodel.Protein);
            if any(tfreg)
                %if strcmp(par_name,'coeff')
                if ~par_valflag
                    if trnmodel.Coefficient(indx,tfreg) ~= 0 ||...
                            trnmodel.RS(indx,tfreg) ~= 0
                        if ~isempty(par_data)
                            if strcmp(par_name,'coeff')%Changing Coefficients                               
                                trnmodel.Coefficient(indx,tfreg) =...
                                                     par_data(ireg);
                                flag = 5;%Coefficient changed to given value
                            elseif strcmp(par_name,'srate')%Changing S Transcription Rates                                
                                trnmodel.srate(indx,tfreg) =...
                                               par_data(ireg);
                            end
                        else
                            if strcmp(par_name,'coeff')%Changing Coefficients                                
                                if trnmodel.RS(indx,tfreg) > 0
                                    trnmodel.Coefficient(indx,tfreg) =...
                                    defparval.accoeff;
                                elseif trnmodel.RS(indx,tfreg) < 0
                                    trnmodel.Coefficient(indx,tfreg) =...
                                    defparval.repcoeff;
                                elseif trnmodel.RS(indx,tfreg) == 2
                                    trnmodel.Coefficient(indx,tfreg) =...
                                    defparval.dualcoeff;
                                end
                                flag = 6;%Coefficient changed to DEFAULT VALUES  
                            elseif strcmp(par_name,'srate')
                                if trnmodel.Coefficient(indx,tfreg) ~= 0 ||...
                                   trnmodel.RS(indx,tfreg) ~= 0
                                    %fprintf('Changing S Transcription Rates\n');
                                    trnmodel.srate(indx,tfreg) =...
                                    defparval.srate;
                                end
                            end
                        end
                    else
                        fprintf('Regulator %s does not regulate Gene %s\n',...
                            trnmodel.Protein{tfreg},...
                            trnmodel.Gene{indx});
                    end
                else
                    if trnmodel.RS(indx,tfreg) > 0
                        trnmodel.Coefficient(indx,tfreg) = defparval.accoeff;
                    elseif trnmodel.RS(indx,tfreg) < 0
                        trnmodel.Coefficient(indx,tfreg) = defparval.repcoeff;
                    elseif trnmodel.RS(indx,tfreg) == 2
                        trnmodel.Coefficient(indx,tfreg) = defparval.dualcoeff;
                    end
                    flag = 6;%Coefficient changed to DEFAULT VALUES
                end
                %end
            else
                fprintf('Regulator %s not found in model\n',regulators{ireg});
            end            
            ireg = ireg + 1;
        end
    end
    
    if strcmp(par_name,'brate')
        trnmodel.brate(indx) = par_data;
%         nreg = 0;
        flag = 1;%Basal Rate changed
    end
    if strcmp(par_name,'drate')
        %drate not yet defined on a gene by gene basis
        %trnmodel.drate(tfg,1) = par_data.drate;
        %So this function will not need to support this parameter
        fprintf('Unsupported parameter\n');
        %fprintf('Only chnage default values\n');
        %defparval.drate = par_data;
        flag = 2; %drate is changed/unchanged
    end
    if strcmp(par_name,'trate')
        %trate is defined only gap. Metabolite associated proteins do not have a trate
        trnmodel.trate(indx,logical(trnmodel.trate(indx,:))) = par_data;
%         nreg = 0;
        flag = 3;%trate is changed
    end
    if strcmp(par_name,'pdrate')
        trnmodel.pdrate(indx) = par_data;
%         nreg = 0;
        flag = 4;%pdrate is changed
    end
    
    %Just use a single coefficient to denote all types of coefficients
    %(activation,inhibition,dual coefficients)
    %coeffs = {'acccoeff';'repcoeff';'dualcoeff'}; 
end

function coeffchange(par_name,indx,par_val,regulators) 
    %Just use a single coefficient to denote all types of coefficients
    %(activation,inhibition,dual coefficients)
    %coeffs = {'acccoeff';'repcoeff';'dualcoeff'};
    
    par_valflag = 0;
    if ~isempty(regulators) && length(regulators) == length(par_val)
        nreg = length(regulators);        
    else%Change parameters for all regulators
        regulators = trnmodel.Regulator(indx);
        nreg = length(trnmodel.Regulator(indx));
        par_valflag = 1;        
    end
    ireg = 1;
    while ireg <= nreg
        tfreg = strcmp(regulators{ireg},trnmodel.Protein);
        if any(tfreg)
            if strcmp(par_name,'coeff')
                if ~par_valflag
                    if trnmodel.Coefficient(indx,tfreg) ~= 0 ||...
                       trnmodel.RS(indx,tfreg) ~= 0
                        trnmodel.Coefficient(indx,tfreg) = par_val(ireg);
                    else
                        fprintf('Regulator %s does not regulate Gene %s\n',...
                                 trnmodel.Protein{tfreg},...
                                 trnmodel.Gene{indx});
                    end
                else
                    if trnmodel.RS(indx,tfreg) > 0
                        trnmodel.Coefficient(indx,tfreg) = defparval.acccoeff;
                    elseif trnmodel.RS(indx,tfreg) < 0
                        trnmodel.Coefficient(indx,tfreg) = defparval.repcoeff;
                    elseif trnmodel.RS(indx,tfreg) == 2
                        trnmodel.Coefficient(indx,tfreg) = defparval.dualcoeff;
                    end
                end                
            end            
        end
        
        ireg = ireg + 1;
    end
end 

% 
% function [fieldindx] = chkfield(Struct,FieldNames)
%         nfields = length(FieldNames);
%         fieldindx = zeros(nfields,1);
%         ifield = 1;
%         while ifield <= nfields
%             fieldindx(ifield) =  isfield(Struct,FieldNames{ifield});               
%             ifield = ifield + 1;
%         end
% end
%     
    


if ~isempty(ProtStruct)
    if isfield(ProtStruct,regulators)
        ProtStruct.allprot = 0;
    else
        ProtStruct.all = 1;
    end
    if isfield(ProtStruct,pname)
        for igene = 1:length(ProtStruct.pname)
        end
    else
        fprintf('No protein name given.');
        fprintf(' Specified parameters will be changed for all genes\n');
    end
    
end




%         pfindx = cellfun(@(x)isfield(x,protfields{ipfield}),...
%                             ProtStruct,'UniformOuput',false);

%             gfindx = cellfun(@(x)isfield(x,genefields{igfield}),...
%                             GeneStruct,'UniformOutput',false);
    



if ~isempty(ProtStruct)
    for iprot = 1:length(ProtStruct)
        pfindx = logical(chkfield(ProtStruct{iprot},ProtFields));
        if ~isempty(ProtFields(pfindx))
            if any(strcmp('trate',ProtFields(pfindx)))
                
                ratechange({'trate'},ProtStruct{iprot});
                
            elseif any(strcmp('pdrate',ProtFields(pfindx)))
                
                ratechange({'pdrate'},ProtStruct{iprot});
                
            end
        end
    end
end        
end
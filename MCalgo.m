%change intial condition as a percentage of input value

%get all input concentrations

%get required metabolite index

%get required percentage shift (+ or -)
    %if specified
        %check if change is applied to specific or all metabolites
            %if specific metabolites
                %add or substract required concentration to or from steady state
            %else
                %change all metabolites
            %end
    %if not specified
        %choose random bounds for concentration changes
            %determine bounds for changes in concentrations    
        %check if change is applied to specific or all metabolites
            %if specific metabolites
                %add or substract required concentration to or from steady state
            %else
                %change all metabolites
            %end 
        
    


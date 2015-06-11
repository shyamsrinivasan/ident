% function [hfig] = savefile(data,datafname,saveData,hfig,figfname)
%**************************************************************************
%Function to save data to text files and figures to png files
%*****Inputs
%data       Data array that is to written to text file. Pass empty array 
%           if no data is to be saved
%datafname  Name of the file in which data is to be written. Existing
%           filenames are overwritten with new data
%saveData   Structure contatining information regarding the directory in 
%           which data is to be stored
%hfig       Handle of figure that is to be saved as a png file. 
%figfname   File name for the png file. No defaults
%*****Outputs
%
%**************************************************************************
%Notes:
%Preceded by writetofile.m - January 2014
function hfig = savefile(data,datafname,saveData,hfig,figfname)
    if nargin < 4
        hfig = 0;
    end   
    if ~iscell(data)
        try
            file_d = [data.t,data.y'];    
        catch
            file_d = [data.t,data.y];
        end
        %New Directory to Store Results
        curr_dir = pwd;
        cd(saveData.dirname);
        if exist(saveData.filename,'dir');
%             fprintf('Directory Already Exists\nWriting Data\n');
        else
            mkdir(saveData.filename);
        end
        cd(saveData.filename);   
        %Write Solution to Text file    
        if ~isempty(datafname)
            fid = fopen(datafname,'w');
            if fid == -1
                fprintf('Cannot Open File\n');
            else
                fclose(fid);
            end
            dlmwrite(datafname,file_d,'delimiter','\t','precision','%5.3e');
        end   
        %Save Figure to Same Directory as above   
        if hfig        
            set(hfig,'PaperPositionMode','auto');
            set(hfig,'InvertHardcopy','off');
            print('-dpng','-r300',figfname);
            set(hfig,'Visible','off');
            %hallfig(1) = hfig;
        else
            hfig = 0;
            %hallfig(1) = 0;
        end
        %whitebg(hsubfig,[0 0 0])
    else
        table = cell2table(data,'VariableNames',...
                {'VariableorParameter','lbRange','ubRange','LB','UB','NewValue'});
        curr_dir = pwd;
        cd(saveData.dirname);
        if exist(saveData.filename,'dir');
            fprintf('Directory Already Exists\nWriting Data\n');
        else
            mkdir(saveData.filename);
        end
        cd(saveData.filename); 
        writetable(table,datafname,'FileType','text','Delimiter','\t');
    end    
    cd(curr_dir);
end
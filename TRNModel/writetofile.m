 %Write Solution to Xcel or Text File
 %Superceded by savefile.m - January 2014
 %
 function status = writetofile(data,newfname)
 if nargin < 2
%      currentfolder = pwd;
%      cd('C:\Users\shyam\Documents\Courses\CHE 1125H\Project\Results');     
%      newfname = datestr(now);
     newfname = sprintf('C:\Users\shyam\Documents\Courses\CHE 1125H\Project\Results\%s.xlsx',...
                datestr(now));
 end
 
 
 %fileid = fopen(newfname,'w');
 
% if fileid == -1
%      fprintf('Results file cannot be opened\n');
%      wflag = fileid;
%      return
% end
 
status = zeros(2,1);
data{2} = data{2}';
 if iscell(data)
%      ntime = length(data{1});
%      [nrows,ncol] = length(data{2});
     status(1) = xlswrite(newfname,data{1},'Sheet1','A1');
     status(2) = xlswrite(newfname,data{2},'Sheet2','A1');
 end
     
%      if ntime == ncol
%          for itime = 1:ntime
%              fprintf('
 %  Solution.initSS{1} = initSS.t;
 %  Solution.initSS{2} = initSS.y;
    
 
 
%  data{1} = Solution.initSS{1} - Time
%  data{2} = Solution.initSS{2} - Expression Values

%fclose(fileid);
 
 end
 
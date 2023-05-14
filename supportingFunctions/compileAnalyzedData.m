function [analyzed_comp, nrSessions] = compileAnalyzedData (fPath, Animal, expID, gratings)

for exp = 1:length(expID)
    if gratings
         dir_path = [fPath filesep Animal filesep 'gratings' filesep expID{exp}];
    else 
         dir_path = [fPath filesep Animal filesep 'barmapping' filesep expID{exp}];
    end
         % get a list of all subfolders
         sessions = dir(dir_path);
         sessions = sessions([sessions(:).isdir]); % keep only directories
         sessions = sessions(3:end); % remove '.' and '..'
         
         % count the number of subfolders whose names contain only numbers
         count = 0;
         for i = 1:length(sessions)
             if regexp(sessions(i).name, '^\d+$')
                 count = count + 1;
             end
         end
         nrSessions(exp) = count;
         
         % compile analyzed mat. file from each session (column) and compile
         % them with different experiments
         for sessNr = 1:nrSessions(exp)
             p = [dir_path filesep num2str(sessNr)];
             cd(p)
             if gratings
                 analyzed_comp(exp, sessNr)= load(['analyzed_passive_' num2str(sessNr) '.mat']);
             else
                 analyzed_comp(exp, sessNr)= load(['analyzed_barmapping_' num2str(sessNr) '.mat'], 'avgStimResponse','avgTempTrace','cont1','roi');
             end
         end
end 
save ([dir_path filesep 'compiled analyzed data.mat'], 'analyzed_comp');


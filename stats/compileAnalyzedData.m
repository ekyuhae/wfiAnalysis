function [analyzed_comp, nrSessions, sessNr_comp] = compileAnalyzedData (fPath, Animal, expID, gratings)
% compile analyzed mat. file from each session and compile them with different experiments
% EK May23

for exp = 1:length(expID)
    if gratings == 1 
         dir_path = [fPath filesep Animal filesep 'gratings' filesep expID{exp}];
    elseif gratings == 0
         dir_path = [fPath filesep Animal filesep 'barmapping' filesep expID{exp}];
    elseif gratings == 2
        dir_path = [fPath filesep Animal filesep 'behavior' filesep expID{exp}]; 
    end
         % get a list of all subfolders
         sessions = dir(dir_path);
         % cfile_comp = "compiled analyzed data.mat";
         % if exist("cfile_comp")
         %     load(cfile_comp)
         sessions = sessions([sessions(:).isdir]); % keep only directories
         sessions = sessions(3:end); % remove '.' and '..'
         s = [];
         count = 0;
         for i = 1:length(sessions)
             if regexp(sessions(i).name, '^\d+$')
             temp = i; % take the session whose name is an integer only
             sessNr = sessions(temp).name;
                 sessNr = str2double(sessNr);

                 p = [dir_path filesep num2str(sessNr)];
                 % if isfolder(p)
                 cd(p)
                 if gratings == 1
                     cfile = ['analyzed_passive_' num2str(sessNr) '.mat'];
                     if exist(cfile)
                         analyzed_comp(exp, sessNr)= load(cfile);
                     end
                 elseif gratings == 0
                     cfile = ['analyzed_barmapping_' num2str(sessNr) '.mat'];
                     if exist(cfile)
                         analyzed_comp(exp, sessNr)= load(cfile, 'avgStimResponse','avgTempTrace','cont1','roi');
                     end
                 elseif gratings ==2
                     cfile = ['analyzed_active_' num2str(sessNr) '.mat'];
                     if exist(cfile)
                         analyzed_comp(exp, sessNr).whole= load(cfile);
                         analyzed_comp(exp, sessNr).hit= load(['analyzed_active_' num2str(sessNr) '_hit.mat']);
                         analyzed_comp(exp, sessNr).miss= load(['analyzed_active_' num2str(sessNr) '_miss.mat']);
                         analyzed_comp(exp, sessNr).CR= load(['analyzed_active_' num2str(sessNr) '_CR.mat']);
                         analyzed_comp(exp, sessNr).LL= load(['analyzed_active_' num2str(sessNr) '_Late.mat']);
                         analyzed_comp(exp, sessNr).FA= load(['analyzed_active_' num2str(sessNr) '_FA.mat']);
                     end
                 end
                 count = count +1;
                 s = [s sessNr];
             end
             sessNr_comp{exp} = s;
         end

         nrSessions(exp) = count;
end 
save ([dir_path filesep 'compiled analyzed data.mat'], 'analyzed_comp');


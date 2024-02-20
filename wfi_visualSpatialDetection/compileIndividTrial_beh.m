function dFFtraces_all = compileIndividTrial_beh (Animal, expID, fPath, n)


    cPath = [fPath Animal filesep 'behavior' filesep expID];
       
         % get a list of all subfolders
         % sessions = dir(cPath);
         % sessions = sessions([sessions(:).isdir]); % keep only directories
         % sessions = sessions(3:end); % remove '.' and '..'
         % 
         % % count the number of subfolders whose names contain only numbers
         % count = 0;
         % for i = 1:length(sessions)
         %     if regexp(sessions(i).name, '^\d+$')
         %         count = count + 1;
         %     end
         % end
         % nrSessions= count;
         % 
         % 
         % compile analyzed mat. file from each session (column) and compile
         % them with different experimentss
         dFFtraces_all = cell(n,2); %binoc and monoc 
         for exp = 1:nrSessions
             p = [cPath filesep num2str(sessNr) filesep];
             cd(p)
             cFile =['compiled analyzed data.mat'];
             if exist(cFile)
                load(cFile);
                
             end 
             dFFtraces_all = cellfun(@(x,y) [x;y], dFFtraces_all, dFFtrace, 'UniformOutput',false);
             dFFtraces_all = cellfun(@(x,y) [x;y], dFFtraces_all, indivTempTrace.binoc, 'UniformOutput',false);
         end 


clear; close all
Animal = 'M230209_1';
hva = [];
% expID_ret = '04-Apr-2023_1'; %'04-Apr-2023_1'; % for retinotopy registration and identifying the visual areas 
 expID_ret = '16-Mar-2023';
lastDay = '12-May-2023_1';


addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
fPath = 'Y:\haider\Data\analyzedData\EKK\WFI';
prev_path = [fPath filesep Animal filesep 'behavior'];
list = dir(prev_path); % changed to dir from ls
explist = []; %list of all dates for the animal (not just new ones given as input)
for i = 1:height(list)
    if double(list(i,1).name(1)) >= 48 && double(list(i,1).name(1) <=57) %makes sure it's only the dates
        explist = [explist;{list(i,1).name}];
    end
end
if ~strcmp(explist{end},lastDay) %cut off list at given last day
     explist(end) = [];
end
cd(prev_path)
load("sessMat.mat")
%%
for d = 1: length(explist)
    expID = explist{d,:};
    %change date format
    date = datestr(expID(1:11), 'yyyy-mm-dd');
    for s = 1:width(sessMat{d})
        sMat = sessMat{d};
        sessNr = sMat(1,s);
        expInfoNr = sMat(2,s);
        wfiBehaviorFunc(Animal, expID, hva, expID_ret, sessNr, num2str(expInfoNr), date, 0, 0, 1)
        fprintf(['finished=' num2str(sessNr) '\n'])
    end
    fprintf(['finished=' expID '\n'])
end
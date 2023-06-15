close all; clear; clc


addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))
addpath(genpath("Y:\haider\Code\behaviorAnalysis"));

Animals = {'M230220_1'}; %cell array
expID = {'26-May-2023'}; %cell array
hva ='RL_AM';
flag = 0 ; %if compiling first time 

for i = 1:length(Animals)
    Animal = Animals{i};
    plotDffTracesAndPeaks_beh_RT('hit', Animal, expID(i), hva,0)
    % plotDffTracesAndPeaks_beh('miss', Animal, expID(i), hva,0)
    plotDffTracesAndPeaks_beh_RT('FA', Animal, expID(i), hva,0)
    %plotDffTracesAndPeaks_beh('CR', Animal, expID(i), hva,0)
    plotDffTracesAndPeaks_beh_RT('LL', Animal, expID(i), hva,0)
end

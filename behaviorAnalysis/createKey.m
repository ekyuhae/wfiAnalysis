clear; close all; clc

animal = 'M230831_1';
keynew = struct;
days = {'2023-10-05' '2023-10-06' '2023-10-07' '2023-10-09'  '2023-10-10'  '2023-10-12'  '2023-10-13' '2023-10-15'  '2023-10-17'...
     '2023-10-19'  '2023-10-20'  '2023-10-23'  '2023-10-25'  '2023-10-26'};
% days = {'2023-10-12'  '2023-10-13' '2023-10-14' '2023-10-15'  '2023-10-17'...
     % '2023-10-19'  '2023-10-20'  '2023-10-23'  '2023-10-25'  '2023-10-26'};

sess =  {[1 3] [1 2 3 4 5 6] [1 2 3 4 5 6 7 8] [2 3 4 5 6 7] [1 2 3 4 5 6] [1 2 3 4 5 6 7] [1 2 3 4 7 8]...
    [1 4 5 6 7 9] [1 2 3 4 5 6 7] [1 3 4 5 7 8 9] [1 2 3 4 7 9 10] [1 2 3 4 5 6] [1 2 3 4 6] [1 2 3 5 6 7]};

for d = 1:length(days)
    keynew.(animal)(d).Date = days{d};
    for s = 1:length(sess{d})
        [~,~,out,beh_all,~, attn] = attention_effects_EK(animal,days{d}, split(num2str(sess{d}))); 
        close all
        keynew.(animal)(d).Sessions(sess{d}(s)).Out = beh_all;
        keynew.(animal)(d).Sessions(sess{d}(s)).Attn = attn;
                keynew.(animal)(d).Out = beh_all;
        keynew.(animal)(d).Attn = attn;

    end
end


%%
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\stats\behavior'];
load("compBeh.mat");
key.(animal) = keynew.(animal);
save("compBeh.mat", "key");
%%

% get number of trials in the compiled data

clear binoc monoc
for exp = 1:size(analyzed_comp,1)
   for s = 1:size(analyzed_comp,2)

        if ~isempty(analyzed_comp(exp,s).hit)
binoc{exp,s} = cellfun(@(x) size(x,1), analyzed_comp(exp,s).CR.indivTempTrace.binoc);
        end 
    end
end

for exp = 1:size(analyzed_comp,1)
    for s = 1:size(analyzed_comp,2)
        if ~isempty(analyzed_comp(exp,s).hit)
monoc{exp,s} = cellfun(@(x) size(x,1), analyzed_comp(exp,s).CR.indivTempTrace.monoc);
        end 
    end
end

%% fg
for exp = 1:size(analyzed_comp,1)
    for s = 1:size(analyzed_comp,2)
        if ~isempty(binoc{exp,s})
            if length(binoc{exp,s}) < 6
                % binoc{exp, s}(end) = NaN;
                % monoc{exp, s}(end) = NaN;
                binoc{exp, s} = [binoc{exp, s}(1) 0 binoc{exp, s}(2:end)];
                monoc{exp, s} = [binoc{exp, s}(1) 0 monoc{exp, s}(2:end)];
            end
        end
    end
end
% c2{1,1} = [0 0 9 7 8 11]
% c2{1,2} = [0 0 10 5 5 10]
%%
sumVectorb = zeros(1,6);
sumVectorm = zeros(1,6);
%%
for i = 1:size(binoc, 1)
   for ii = 1:size(binoc, 2)
    if ~isempty(binoc{i, ii})
        sumVectorb = sumVectorb + binoc{i, ii};
    end
    end
end

for i = 1:size(monoc, 1)
   for ii = 1:size(monoc, 2)
    if ~isempty(monoc{i, ii})
        sumVectorm = sumVectorm + monoc{i, ii};
    end
    end
end
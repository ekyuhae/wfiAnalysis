function dff = normDFF (data, baseline)
% EK Jan23
% compute df/f 
data = single(data); 

if numel(size(data)) == 3
    [A,B,C] = size(data);
    dataAvg = mean(data(:,:,baseline),3);
    data = bsxfun(@minus, data, dataAvg); % subtract baseline mean
    data = bsxfun(@rdivide, data, dataAvg); % divide by baseline mean
    dff = reshape(data,[],C);
else
    dataAvg = mean(data, baseline);
    data = bsxfun(@minus, data, dataAvg); % subtract baseline mean
    dff = bsxfun(@rdivide, data, dataAvg); % divide by baseline mean
end 
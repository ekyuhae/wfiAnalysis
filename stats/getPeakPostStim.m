function [peak_value, peak_index] = getPeakPostStim(data, framerate, window_duration, start_time, peakOnly)
% data = temporal trace vector (row vector)
% framerate = effective frame rate (e.g. 15Hz for blue channel data)
% window_duration = post stimulus window in second e.g. 1
% start_time = onset time in second e.g. 2

% Compute the number of samples in the time window
window_samples = floor(window_duration * framerate);

% Compute the starting and ending indices of the time window
start_index = round(start_time * framerate)+1;
end_index = start_index + window_samples - 1;

% Check if the window extends beyond the end of the data vector
if end_index > length(data)
    error('Window extends beyond the end of the data vector');
end

% Find the maximum value within the specified time window
[peak_value, peak_index] = max(data(start_index:end_index));

if isempty(peakOnly)
    peak_value = peak_value - data(start_index-1);  % calculate the evoked response's magnitude
else
    peak_value = peak_value;
end

% Add the starting index of the window to get the peak index relative to the whole data vector
peak_index = peak_index + start_index - 1;

end

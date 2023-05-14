function [peak_value, peak_index] = getPeakPostStim(data, framerate, window_duration, start_time)

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
peak_value = peak_value - data(start_index-1);

% Add the starting index of the window to get the peak index relative to the whole data vector
peak_index = peak_index + start_index - 1;

end

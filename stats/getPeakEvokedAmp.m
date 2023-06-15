function magnitude = getPeakEvokedAmp(RT, ST, framerate, stimOnset)

% RT: lickTrigged 3d temporal trace (N_contrast x N_timepoints x
% N_sessions)
% ST: stimTrigged 3d temporal trace (N_contrast x N_timepoints x
% N_sessions)
% framerate = effective frame rate (e.g. 15Hz for blue channel data)
% window_duration = post stimulus window in second e.g. 1
% start_time = onset time in second e.g. 2

start_time = round(stimOnset * framerate)+1;

for iSess = 1: size(RT,3)
    for iCont = 1: size (RT,1)
        % [peakPreLick, idx] = getPeakPostStim (RT(iCont,:,iSess), framerate, preLickWindow, 0, 1);
        [peakPreLick, peak_index] = max(RT(iCont,start_time-1:-1:start_time-10 ,iSess)); % lick within .7s window, so find peak within 10 frames
        baseline = ST(iCont,start_time, iSess);
        magnitude(iCont,iSess)= peakPreLick - baseline;
    end
end

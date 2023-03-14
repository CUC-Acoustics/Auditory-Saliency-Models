function pitch = pitch_feature(in,fs)
% This function is used to compute the one-dimensional pitch feature of the
% audio.
% INPUT
% -- in: Audio signal vector.Audio signal vector.
% -- fs: Sampling frequency.

audio = in;
Fs = fs;
audio = audio./max(abs(audio));

% Design a bandpass filter to filter out noise interference.
Fpass1 = 60;   % Passband starting frequency.
Fpass2 = 2000;  % Passband cutoff frequency.
filterOrder = 3;  % Filter order.
[b, a] = butter(filterOrder, [Fpass1, Fpass2]/(Fs/2), 'bandpass');
audio = filter(b, a, audio);

% Frame segmentation parameters.
frame_size = round(0.02 * Fs); 
frame_shift = round(0.02 * Fs);
min_pitch = 60;
max_pitch = 2000;

% Compute pitch
[pitch, ~] = calculate_pitch(audio, frame_size, frame_shift, min_pitch, max_pitch, Fs);


function [pitch, frames] = calculate_pitch(audio, frame_size, frame_shift, min_pitch, max_pitch, fs)
% INPUT:
%   - audio: Audio signal vector.
%   - frame_size: Length of each frame (in number of samples).
%   - frame_shift: Overlap length between frames (in number of samples).
%   - min_pitch: Minimum pitch (in Hz).
%   - max_pitch: Maximum pitch (in Hz).
% OUTPUT:
%   - pitch: Pitch vector of each frame (in Hz).
%   - frames: Frame num

% Segment the audio signal into frames.
num_samples = length(audio);
num_frames = floor((num_samples - frame_size) / frame_shift) + 1;
frames = zeros(frame_size, num_frames);
for i = 1:num_frames
    start = (i-1)*frame_shift + 1;
    frames(:,i) = audio(start:start+frame_size-1);
end

% Set an energy threshold to avoid noise interference.
energy = zeros(1, num_frames);
for i = 1:num_frames
    energy(i) = sum(frames(:,i).^2);
end
energy = energy./max(energy);


% Compute the pitch for each frame.
pitch = zeros(1, num_frames);
for i = 1:num_frames
    frame = frames(:,i);
    r = xcorr(frame);
    r = r(frame_size:end);
    max_lag = ceil(fs/min_pitch);
    min_lag = ceil(fs/max_pitch);
    r = r(min_lag:max_lag);
    [pks,locs] = findpeaks(r);
    if isempty(pks) || energy(i) < 0.1
        pitch(i) = 0;
    else
        [~,index] = max(pks);
        pitch(i) = 1 / (min_lag + locs(index) - 1);
    end
end
end


end

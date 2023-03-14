function out = waveform(in,fs)
% This function is used to compute the envelope feature of the input audio.
%
% INPUT
% -- in: Audio signal vector.Audio signal vector.
% -- fs: Sampling frequency.

h1 = abs(hilbert(in(:,1)));
[b,a]=butter(6,60/(fs/2),'high');

y = abs(filter(b,a,h1));
r        = 10;
sigma    = 1;
fram=fs*0.02;
yG = Gaussianfilter(r, sigma, y);
LyG = ceil(length(yG)/fram);
out = zeros(LyG*fram,1);
for i=1:length(yG)
    out(i)=yG(i);
end
out=enframe(out,fram,fram);
out=mean(out,2)';
out=out/max(out);
return


function y_filted = Gaussianfilter(r, sigma, y)

% Generate a one-dimensional Gaussian filter template.
GaussTemp = ones(1,r*2-1);
for i=1 : r*2-1
    GaussTemp(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end

% Gaussian filtering.
y_filted = conv(y,GaussTemp);
return


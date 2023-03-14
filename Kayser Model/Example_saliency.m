

% scipt for simple testing of auditory saliency map
% see figure S1 in the paper
%
% create a short and a long tone on background noise

FS =16000; % sampling rata
Dur = 2.2; % total duration sec

x = [1/FS:1/FS:Dur];
Bgnd = randn(1,length(x))*0.5; % background sound - noise

% short tone
Freq0 = 2000; % Hz.
% make a hanning windowed mask
d = round(10*FS/1000);
h = hanning(d);
T0 = FS*0.5;
X = zeros(size(x));
X(T0:T0+30*FS/1000) = 1;
X(T0:T0+d/2-1) = h(1:d/2);
X(T0+30*FS/1000-d/2+1:T0+30*FS/1000) = h(d/2+1:end);
% create the tone and multiply by mask
Tone1 = sin(2*pi*x*Freq0).*X; % first tone

% longer tone
d = round(40*FS/1000);
h = hanning(d);
X = zeros(size(x));
T0 = FS+FS/4;
X(T0:T0+500*FS/1000) = 1;
X(T0:T0+d/2-1) = h(1:d/2);
X(T0+500*FS/1000-d/2+1:T0+500*FS/1000) = h(d/2+1:end);
% create the tone and multiply by mask
Tone2 = sin(2*pi*X.*x*Freq0).*X;  % second tone

Sound1 = Tone1+Tone2+Bgnd*1.2;  %put tones and backgrpound together

[S1,f,t1] = specgram(Sound1,1024,FS,800,778); % compute spectrogram of this sound

S1 =log(abs(S1)); % make intensity map

SALIENCY = Saliency_map(S1,4); % compute saliency map

SAL = SALIENCY.eo + SALIENCY.esi + SALIENCY.epi; % combine saliency maps from the three different filters

figure;
imagesc(SAL);
title('Saliency map');


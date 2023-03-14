
clearvars

% Reading audio
[s,fs] = audioread('audio_test.wav');
s = s./max(abs(s));

% Cochlear filtering 
%
% In the "scm" input argument, the parameters are as follows
% 128: the number of channels 
% 100: the lowest frequency
% 8000: the highest frequency 
% "5": time resolution is 5 milliseconds 
% "0": without high pass filter in cochlear output
[eResp,fx,cf,tx] = scm(s(:,1),fs,[128 100 8000],5, 0);

% Lateral inhibition network.
S1 = LIN(eResp,1);

% Compute saliency maps with Gaussian pyramid sampling at 4 levels.
SALIENCY = Kalinli_Saliency_map(S1,4);

% Linearly combine the feature maps.
SAL = SALIENCY.eo + SALIENCY.esi + SALIENCY.epi + SALIENCY.o1 + SALIENCY.o2 + SALIENCY.p;

% Draw the saliency map.
figure(1);
ndx = linspace(0,(length(s)/fs)*1000,size(SAL,2));
mdx = linspace(0,floor(fs/2), size(SAL,1));
imagesc(ndx,mdx,SAL);
title('Saliency map');
xlabel('Time/ms')
ylabel('Frequency/Hz')



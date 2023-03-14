
clearvars

% Reading audio
[s,fs] = audioread('audio_test.wav');

% calculate envelope features
KA.data{1}{1} = waveform(s,fs);

% calculate pitch features
KA.data{2}{1} = pitch_feature(s,fs);

% calculate frequency features
KA.data{3}{1} = frequency(s);

% calculate bandwidth and rate features
[KA.data{4}{1},KA.data{5}{1}] = rs_feature(s);

% resampling to obtain different scales
% "4" is the scale of  center-surround difference

N=4;
for i=1:5
    for j=1:N
        KA.data{i}{j} = imresize(KA.data{i}{1}, 1/(2^(j-1)));
    end
end
for i=1:5
    Mapsize=size(KA.data{i}{1})/4;
    for j=1:N
        KA.data{i}{j} = imresize(KA.data{i}{j},Mapsize);
    end
end


% the range of  center-surround difference
scale_interact =[1,2];
ONOFF = [1,0];

% center-surround difference
KA = CenterSurroundPyramid(KA,scale_interact,N,ONOFF);

% Perform normalization processing.
SALIENCY.wave = zeros(size(KA.data{1}{1}));
SALIENCY.p = zeros(size(KA.data{2}{1}));
SALIENCY.freq = zeros(size(KA.data{3}{1}));
SALIENCY.r = zeros(size(KA.data{4}{1}));
SALIENCY.s= zeros(size(KA.data{5}{1}));



% envelope
for n=1:size(KA.data{1},2)
    if max(KA.data{1}{n}) ~=0
          SALIENCY.wave = SALIENCY.wave + KA.data{1}{n}/max(KA.data{1}{n});
    else
      SALIENCY.wave = SALIENCY.wave + KA.data{1}{n};
    end
end
  
% pitch
for n=1:size(KA.data{2},2)
    if max(KA.data{2}{n}) ~=0
          SALIENCY.p = SALIENCY.p + KA.data{2}{n}/max(KA.data{2}{n});
    else
      SALIENCY.p = SALIENCY.p + KA.data{2}{n};
    end
end
 
% frequency
for n=1:size(KA.data{3},2)
    if max(KA.data{3}{n}) ~=0
          SALIENCY.freq = SALIENCY.freq +DOG(KA.data{3}{n},20);
    else
          SALIENCY.freq = SALIENCY.freq + KA.data{3}{n};
    end
end
  
% rate
for n=1:size(KA.data{4},2)
    if max(KA.data{4}{n}) ~=0
          SALIENCY.r = SALIENCY.r +DOG(KA.data{4}{n},20);
    else
          SALIENCY.r = SALIENCY.r + KA.data{4}{n};
    end
end
  
% scale
for n=1:size(KA.data{5},2)
    if max(KA.data{5}{n}) ~=0
          SALIENCY.s = SALIENCY.s +DOG(KA.data{5}{n},20);
    else
          SALIENCY.s = SALIENCY.s + KA.data{5}{n};
    end
end

% Calculate the mean of the two-dimensional features along the frequency
% dimension to convert them into one-dimensional features.
SALIENCY.freq = mean(SALIENCY.freq,1);
SALIENCY.r = mean(SALIENCY.r,1);
SALIENCY.s = mean(SALIENCY.s,1);

% Set the maximum value of all one-dimensional features to 1.
if max(SALIENCY.wave)~=0
    SALIENCY.wave = SALIENCY.wave/max(SALIENCY.wave);
else
    SALIENCY.wave = SALIENCY.wave;
end

if max(SALIENCY.p)~=0
    SALIENCY.p = SALIENCY.p/max(SALIENCY.p);
else
    SALIENCY.p = SALIENCY.p;
end

SALIENCY.freq = SALIENCY.freq/max(SALIENCY.freq);
SALIENCY.r = SALIENCY.r/max(SALIENCY.r);
SALIENCY.s = SALIENCY.s/max(SALIENCY.s);

% Linearly combine all features.
Saliency = SALIENCY.wave + SALIENCY.p + SALIENCY.freq + SALIENCY.r + SALIENCY.s;
ndx = linspace(0,(length(s)/fs)*1000,numel(Saliency));
plot(ndx,Saliency);
xlabel('Time/ms')
title('Saliency map')





%%%%%%%%%%%%%%%%%%%%%%%%
function CS = CenterSurroundPyramid(PYR,offsets,N,ONOFF)
% center surround interactions on different spatial scales
% offsets: indicates the differences in spatial scale used for the
%          interaction
% for each feature map
cnt=1;
for n1=1:3
  for n2=n1+offsets
    if n2<=N
      % loop the features
      for f=1:length(PYR.data)
        if ONOFF(1)
          % on-off
          dummy = PYR.data{f}{n2}-PYR.data{f}{n1};
          CS.data{f}{cnt} =  dummy.*(dummy>0);
         CS.helper{f}(cnt) = abs(n2-n1);
        end
        if ONOFF(2)
          % off-on
          dummy = PYR.data{f}{n2}-PYR.data{f}{n1};
          CS.data{f}{ZZ2,cnt} =  dummy.*(dummy>0);
          CS.helper{f}(cnt) = abs(n2-n1);
        end
      end
      cnt=cnt+1;
    end
  end
end

return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function M = DOG(in,time)
% Perform normalization using a Dog filter.

in = in/max(max(in));
in=in.*(in>0);
[~,yw]=size(in);
x=0;
for ii = linspace(-20,20,41)
    x = x + 1;
    y = 0;
    for jj = linspace(-20,20,41)
        y = y + 1;
        fex=0.02*yw;
        fin=0.25*yw;
        cex=0.5;
        cin=1.5;
        dog(x,y)=-(((cex^2/2*pi*fex^2)*exp(-(ii^2+jj^2)/2*fex^2))-((cin^2/2*pi*fin^2)*exp(-(ii^2+jj^2)/2*fin^2)));
    end
end

M=in;
for r=1:time
    M= M + conv2(M,dog,'same')-0.02;
    M=M.*(M>0);
end
return
end

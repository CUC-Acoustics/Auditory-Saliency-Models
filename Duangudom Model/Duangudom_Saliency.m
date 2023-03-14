
clearvars

% Reading audio
[x,fs]=audioread('audio_test.wav');

% The functions "wav2uad" and "aud2cor" are both from the NLS toolbox,
% which is provided in this folder. "wav2uad" utilizes the early auditory
% processing system of the cochlear model to obtain the auditory
% spectrogram.
%
% In the "paras" input argument, the parameters are as follows
% 20: frame length is 20ms
% 128: time constant is 128ms
% -2: no nonlinear processing is performed
% 0: sampling rate of 16kHz
paras = [20,128,-2,0];
y = wav2aud(x, paras,'p_o');

% rv: the frequency range of the time filter,
% sv: represents the scale range of the frequency filter.
rv=linspace(0.5,32,14);  %rate vector in Hz
sv=linspace(0.25,8,6); %scale vector in cyc/oct

% "aud2cor"  processes the auditory spectrogram through the A1 cortical
% layer and can store the output as a 'cr.cor' file or cr. The variable
% size(cr) is a 4-dimensional array. The first two dimensions represent the
% size of the frequency scale (sv) and twice the size of the time scale
% (2*rv). The last two dimensions represent the size of the auditory
% spectrogram.
cr = aud2cor(y, paras, rv, sv, 'cr.cor');
[k1,k2,k3,k4]=size(cr);
k2=k2/2;
for i=1:k1
    for j=1:k2
        cr0=cr(i,j,:,:);
        cr0=squeeze(cr0)';
        SL{i}{j}=abs(cr0);
        SL{i}{j}=normalizemap(SL{i}{j},50);
    end
end


% Global Energy
SALIENCY.p = zeros(size(SL{1}{1}));
for i=1:3
    for j=1:8
        Du{1}{i}{j}=SL{i}{j};
        SALIENCY.p=SALIENCY.p+Du{1}{i}{j};
    end
end

% Temporal Modulation
SALIENCY.t = zeros(size(SL{1}{1}));
sti=0;
for i=1:3
    sti=sti+1;
    stj=0;
    for j=9:14
        stj=stj+1;
        Du{2}{sti}{stj}=SL{i}{j};
        SALIENCY.t =SALIENCY.t +Du{2}{sti}{stj};
    end
end
% Spectral Modulation
SALIENCY.f = zeros(size(SL{1}{1}));
sfi=0;
for i=4:6
    sfi=sti+1;
    sfj=0;
    for j=1:8
        sfj=sfj+1;
        Du{3}{sfi}{sfj}=SL{i}{j};
        SALIENCY.f=SALIENCY.f+Du{3}{sfi}{sfj};
    end
end

% High Temporal & Spectral Modulation
SALIENCY.ft = zeros(size(SL{1}{1}));
stfi=0;
for i=4:6
    stfi=stfi+1;
    stfj=0;
    for j=9:14
        stfj=stfi+1;
        Du{4}{stfi}{stfj}=SL{i}{j};
        SALIENCY.ft=SALIENCY.ft+Du{4}{stfi}{stfj};
    end
end

% Normalization
SALIENCY.p = normalizemap(SALIENCY.p,50);
SALIENCY.t = normalizemap(SALIENCY.t,50);
SALIENCY.f = normalizemap(SALIENCY.f,50);
SALIENCY.ft = normalizemap(SALIENCY.ft,50);

Saliency=SALIENCY.p+SALIENCY.t+SALIENCY.f+SALIENCY.ft;

% Draw the saliency map.
[M,N]=size(Saliency);
ndx = (1:N) * paras(1);
mdx = (1:M) * fs / 2 / M;
imagesc(ndx,mdx,Saliency);
title('Saliency map');
xlabel('Time/ms');
ylabel('Frequency/Hz')



%%%%%%%%%%%%%%%%%
% Normalization
function out = normalizemap(in,win)
% function that normalizes a feature map with respect to 
% local maxima. In the visual case, this is a glocal spatial operation. 
% Here we need an interaction global in frequency but more localized in temporal
% dimension.
% use a sliding window analysis that normalizes windows of length win
% independently, but taking the data from a 3*win window into account.

warning off 
ntp = size(in,2);
winborder = [1:win:ntp];
out = zeros(size(in));
if mod(ntp,win)~=0
  winborder = [winborder,ntp];
end

% before normalizing, apply a mask to avoid edge-effects

scale = 11;
h = hanning(scale*2);
mask1 = ones(size(in));
mask2 = mask1;
mask1(1:scale,:) = h(1:scale)*ones(1,size(mask1,2));
mask1(end-scale+1:end,:) = h(scale+1:end)*ones(1,size(mask1,2));
mask2(:,1:scale) = ones(size(mask1,1),1)*(h(1:scale)');
mask2(:,end-scale+1:end) = ones(size(mask1,1),1)*(h(scale+1:end)');

in = in.*mask1.*mask2;
% get the local extrema
in = in-min(in(:));
in = in./max(in(:));
[LocMa,LocMi] = localextrema(in);

for W=1:length(winborder)-1
  I_win = [winborder(W):winborder(W+1)];
  I_all = round([winborder(W)-win:winborder(W+1)+win/3]);
  % clip this interval to the data range
  I_all = I_all(find(I_all>0));
  I_all = I_all(find(I_all<=(ntp-2)));
  data = in(:,I_all);
  max_data = max(data(:));
  data = data./max_data;
  globmax = max(data(:));
  
  LocMaX = find(LocMa(:,I_all)); 
  LocMaX = data(LocMaX(:));
  % cancel the global maxima from this list.
  LocMaX = LocMaX(find((LocMaX~=1)));
  LocMaX = LocMaX./max_data;
  % now normalize the data in the smaller interval
  if isempty(LocMaX)
    % no local maximum left. Don't normalize but report
    fprintf('problem with normalization\n');
%     keyboard;
    out_n = (in(:,I_win))./max_data;
  else
    out_n = (in(:,I_win))./max_data;
    out_n = out_n*((1-mean(LocMaX))^2);
  end
  out(:,I_win) = out_n;
end
warning on;
return;
end

function [Maxima,Minima] = localextrema(in)

% function [Maxima,Minima] = localextrema(in)
% 
% find the local maxima and minima of a function

d1 = diff(in,1,1);
d2 = diff(in,1,2);
% find zeros by detecting zero crossings
dum1 = d1(1:end-1,:);
dum2 = d1(2:end,:);
% downwards
cross11 = (dum1>0).*(dum2<=0);
% upward crossing
cross21 =(dum1<0).*(dum2>=0);
dum1 = d2(:,1:end-1);
dum2 = d2(:,2:end);
% downwards
cross12 =(dum1>0).*(dum2<=0);
% upward crossing
cross22 = (dum1<0).*(dum2>=0);

s = size(in)-2;
cross11 = cross11([1:s(1)],[1:s(2)]);
cross12 = cross12([1:s(1)],[1:s(2)]);
cross21 = cross21([1:s(1)],[1:s(2)]);
cross22 = cross22([1:s(1)],[1:s(2)]);

cross11 = [zeros(s(1),1),cross11];
cross12 = [zeros(s(1),1),cross12]; 
cross11 = [zeros(1,s(2)+1);cross11];
cross12 = [zeros(1,s(2)+1);cross12];

cross21 = [zeros(s(1),1),cross21];
cross22 = [zeros(s(1),1),cross22];
cross21 = [zeros(1,s(2)+1);cross21];
cross22 = [zeros(1,s(2)+1);cross22];

% a local maximum occurs if both derivatives cross zero downwards
Maxima = cross11.*cross12;
Minima = cross21.*cross22;
return

end

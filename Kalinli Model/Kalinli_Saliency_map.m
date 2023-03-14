function SALIENCY = Kalinli_Saliency_map(img,N)

% the range of scales at which center surround interactions are computed
scale_interact =[2,3];
% Use on_off and off_on maps
ONOFF = [1,0];

% prepare images ------------------------------------------------------------------

% insert an image border to avoid edge artifacts
img2 = replicate_img(img);

% resample spectrogram on different scales 
for n=1:N
  IMG{n} = imresize(img2,1/(2^(n-1)),'nearest');%imresize(A,m)
end

% the size at which the maps are stored 
MapSize = ceil(size(img)/2);

% extract the features on first order --------------------------------------
LEV = 1;  % the level of the pyramid. 1 - filter responses

for n=1:N
    
  % Compute intensity feature (EO)
  Filter = compute(pi ,0);
  PYR.data{1}{LEV}{n} = abs(filter2(Filter,IMG{n},'same')); % Convolve IMG{n} with the filter, and return a result with the same size as IMG{n}.
  PYR.data{1}{LEV}{n} = imresize(cropimg(PYR.data{1}{LEV}{n}),MapSize,'nearest'); % After cropping the border, resize to a size of MapSize.

  
  % Compute frequency contrast feature (ESI)
  Filter = compute(pi/2 ,2);
  PYR.data{2}{LEV}{n} = abs(filter2(Filter,IMG{n},'same'));
  PYR.data{2}{LEV}{n} = imresize(cropimg(PYR.data{2}{LEV}{n}),MapSize,'nearest');
  
  % Compute Temporal contrast feature (EPI)
  Filter = compute(pi ,1);
  PYR.data{3}{LEV}{n} = abs(filter2(Filter,IMG{n},'same'));
  PYR.data{3}{LEV}{n} = imresize(cropimg(PYR.data{3}{LEV}{n}),MapSize,'nearest');
  
  % Compute orientations feature £¨0£©
  Filter = compute(4/pi ,2);
  PYR.data{4}{LEV}{n} = abs(filter2(Filter,IMG{n},'same'));
  PYR.data{4}{LEV}{n} = imresize(cropimg(PYR.data{4}{LEV}{n}),MapSize,'nearest');
  
  Filter = compute(3*pi/4 ,2);
  PYR.data{5}{LEV}{n} = abs(filter2(Filter,IMG{n},'same'));
  PYR.data{5}{LEV}{n} = imresize(cropimg(PYR.data{5}{LEV}{n}),MapSize,'nearest');
  
  % Compute pitch distribution feature £¨P£©
  PYR.data{6}{LEV}{n} = er_pitch(IMG{n} , img);
  PYR.data{6}{LEV}{n} = imresize(cropimg(PYR.data{6}{LEV}{n}),MapSize,'nearest');
  
  % report if pyramid is too high (i.e. if filters get to big)
  if sum(size(Filter)>size(IMG{n}/2))
    fprintf('Filter size has reached half the data size\n');
    fprintf('%d and %d',size(Filter,1),size(IMG{n},1)/2);
  end
end

% center surround interactions on different spatial scales ---------------
PYR = CenterSurroundPyramid(LEV,PYR,scale_interact,N,ONOFF);

% Linearly combine the feature maps of different scales and normalize using
% a Difference of Gaussians (DoG) filter.
LEV = 2;
SALIENCY.eo = zeros(size(PYR.data{1}{LEV}{1}));
SALIENCY.esi = zeros(size(PYR.data{2}{LEV}{1}));
SALIENCY.epi = zeros(size(PYR.data{3}{LEV}{1}));
SALIENCY.o1 = zeros(size(PYR.data{4}{LEV}{1}));
SALIENCY.o2= zeros(size(PYR.data{5}{LEV}{1}));
SALIENCY.p= zeros(size(PYR.data{6}{LEV}{1}));

for O=find(ONOFF)
  % EO
  for n=1:size(PYR.data{1}{LEV},2) 
    SALIENCY.eo = SALIENCY.eo + DOG(PYR.data{1}{LEV}{O,n},20);
  end
  
  % ESI
  for n=1:size(PYR.data{2}{LEV},2)
    SALIENCY.esi = SALIENCY.esi + DOG(PYR.data{2}{LEV}{O,n},20);
  end
  
  % EPI)
  for n=1:size(PYR.data{3}{LEV},2)
    SALIENCY.epi = SALIENCY.epi + DOG(PYR.data{3}{LEV}{O,n},20);
  end
  
  % o1)
  for n=1:size(PYR.data{4}{LEV},2)
    SALIENCY.o1 = SALIENCY.o1 + DOG(PYR.data{4}{LEV}{O,n},20);
  end
  
  %o2
  for n=1:size(PYR.data{5}{LEV},2)
    SALIENCY.o2 = SALIENCY.o2 +DOG(PYR.data{5}{LEV}{O,n},20);
  end
  
  %p
  for n=1:size(PYR.data{6}{LEV},2)
    SALIENCY.p = SALIENCY.p +DOG(PYR.data{6}{LEV}{O,n},20);
  end
  
end

SALIENCY.eo = SALIENCY.eo / (n*length(find(ONOFF)));
SALIENCY.esi= SALIENCY.esi / (n*length(find(ONOFF)));
SALIENCY.epi = SALIENCY.epi / (n*length(find(ONOFF)));
SALIENCY.o1 = SALIENCY.o1 / (n*length(find(ONOFF)));
SALIENCY.o2 = SALIENCY.o2 / (n*length(find(ONOFF)));
SALIENCY.p = SALIENCY.p / (n*length(find(ONOFF)));

return;

%%%%%%%%%%%%%Local function%%%%%%%%%%%


function gabor_k = compute(theta ,m)
t=3;
g=0.56*t;
r=g;
x = 0;
f0 = 0.2; 
if m==0
    g = 0.5;
    r = 0.3;
end
for i = linspace(-8,8,11)
    x = x + 1;
    y = 0;
    for j = linspace(-8,8,11)
        y = y + 1;
        x1 = i*cos(theta) + j*sin(theta);
        y1 = -i*sin(theta) + j*cos(theta);
        gabor_k (y,x)= exp(-(f0^2*x1^2/(2*r^2)+f0^2*y1^2/(2*g^2)))*cos((2*pi*f0*x1)/t); 
    end
end
if m ==1
    [~,yw]=size(gabor_k);
    gabor_k(:,1:ceil(yw/3))=0;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PYR = CenterSurroundPyramid(LEVin,PYR,offsets,N,ONOFF)
% center surround interactions on different spatial scales
% offsets: indicates the differences in spatial scale used for the
%          interaction
LEV = LEVin+1;
% for each feature map
cnt=1;
for n1=1:N
  for n2=n1+offsets
    if n2<=N
      % loop the features
      for f=1:length(PYR.data)
        if ONOFF(1)
          % on-off
          dummy = PYR.data{f}{LEVin}{n1}-PYR.data{f}{LEVin}{n2};
          PYR.data{f}{LEV}{1,cnt} =  dummy.*(dummy>0);
          PYR.helper{f}{LEV}(1,cnt) = abs(n2-n1);
        end
        if ONOFF(2)
          % off-on
          dummy = PYR.data{f}{LEVin}{n2}-PYR.data{f}{LEVin}{n1};
          PYR.data{f}{LEV}{ZZ2,cnt} =  dummy.*(dummy>0);
          PYR.helper{f}{LEV}(cnt) = abs(n2-n1);
        end
      end
      cnt=cnt+1;
    end
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = er_pitch(in, orig_in)
% This function is used to extract pitch features from a two-dimensional
% auditory spectrogram. For more details, please refer to the paper "On the
% importance of time a temporal representation of sound" in this folder.

% Here, "in" refers to the input auditory spectrograms at different scales,
% and "or_in" refers to the original size of the auditory spectrogram. This
% is used to ensure that the different scale feature maps have different
% frame lengths.
[xw,yw] = size(in);
[~,or_yw] = size(orig_in);

% Because the time resolution was set to 5ms earlier, "8" here represents a
% duration of 40ms. This means that short-term autocorrelation is
% calculated over a length of 40ms.
n = ceil((8)/(or_yw/yw));
Mflag = floor(n/2);
p = zeros(xw,floor(yw/n)*(2*Mflag+1));

for j = 1:xw
    t = in(j,:);
    for m = 1:floor(yw/n)
        xr = xcorr(t((m-1)*n+1 : m*n),Mflag);
        p(j,(m-1)*(2*Mflag+1)+1 : m*(2*Mflag+1)) = xr;
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = cropimg(img)
% crop the image border. Undo replicate_img
% the input consists of 0.5,1,0.5 times the real image.
s = ceil(size(img)/4);
out = img(s(1)+1:end-s(1)+1,s(2)+1:end-s(2)+1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = replicate_img(img)
% replicate the image to its borders to later avoid
% edge artifacts. This is undone by cropimg
out = zeros(size(img,1)*3,size(img,2)*3,size(img,3));
S = size(img);
if length(S)==2
  S(3) = 1;
end

for s=1:S(3)
  % top
  k=1; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=3; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=2; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = (flipud(img(:,:,s)));
  % bottom
  k=2;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = flipud(img(:,:,s));
  k=1;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=3;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  % left right
  k=1;  out([1:S(1)]+(S(1)),[1:S(2)],s) = fliplr(img(:,:,s));
  k=3;  out([1:S(1)]+(S(1)),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(img(:,:,s));
  out([1:S(1)]+S(1),[1:S(2)]+S(2),s) = img(:,:,s);
end
out = out(ceil(S(1)/2):end-ceil(S(1)/2)+1,ceil(S(2)/2):end-ceil(S(2)/2)+1,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = DOG(in,time)
% This is the normalization function. First rescale the values to be within
% the range of 0 to 1, then set any negative values to 0.
in = in/max(max(in));
in=in.*(in>0);
[~,yy]=size(in);
x=0;
for i = linspace(-10,10,21)
    x = x + 1;
    y = 0;
    for j = linspace(-10,10,21)
        y = y + 1;
        fex=0.02*yy;
        fin=0.25*yy;
        cex=0.5;
        cin=1.5;
        dog(x,y)=-(((cex^2/2*pi*fex^2)*exp(-(i^2+j^2)/2*fex^2))-((cin^2/2*pi*fin^2)*exp(-(i^2+j^2)/2*fin^2)));
    end
end
% Perform Gaussian difference.
M=in;
for r=1:time
    M= M + conv2(M,dog,'same')-0.02;
    M=M.*(M>0);
end
return

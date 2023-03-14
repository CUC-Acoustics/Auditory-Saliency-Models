

% load Sound stimuli 
load('Sounds_saliencyexp_1.mat','Sound')
% Sound is a 52x3 cell array with the sound samples, consisting of 52 sounds
% presented on 3 different backgrounds, i.e. Sound{sample,background}.
% During the experiments pairs of randomly chosen sounds were presented and subjects
% had to indicate which one was more salient. In between catch trials were presented,
% where the same sound was presented twice.


% load Saliency ratings from model and subjects
load('Data_saliencyexp_exp1.mat','TOTAL');
% Total is a 7 cell array, containing the data for each of the 7 subjects.
% 1st Dimension: Individuals Trials (different n for each subject)
% the 2nd dimension contains: 
%    colums 1-3 : [index of stimulus 1 (i.e. 1-52), index of stimulus 2 (i.e. 1-52), the background (i.e. 1-3)]
%    column 4: the rating of the subject, with 
%    column 5: the saliency difference obtained from the saliency map (this is the numerical differnece between the
%                      numerical saliency values for stimulus 2 - stimulus 1)
%    column 6: the 'saliency' difference obtained from the sound intensity (spectrogram) (this is the numerical 
%                      differnece between the numerical saliency values for stimulus 2 - stimulus 1)


% How to make sense of this data: The subjects indicated which sound is more salient (response 1 or 2), or indicated
% that both sounds are equally salient (response 0). If sound 2 is more salient than sound 1, according to the model,
% then the numerical saliency differnece (column 5) is positive. An agreement of subjects rating and model
% hence predicts that when column 5 is positive, the subjects response should be 2. In other words, the correlation 
% between columns 4 and 5 should be high.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis

% compute correlations between subjects rating and model

for S=1:7
   % find non-catch trials:
   ind = find(TOTAL{S}(:,4)~=0);
   rating = TOTAL{S}(ind,4);
   model = TOTAL{S}(ind,5);
   % make model prediction binary:
   model = double(model>0);
   cc = corrcoef(model,rating);
   Corr_saliency(S) = cc(1,2)
   
end
% do the same for the predictions based on sound intensity
for S=1:7
   % find non-catch trials:
   ind = find(TOTAL{S}(:,4)~=0);   
   rating = TOTAL{S}(ind,4);
   model = TOTAL{S}(ind,6);
   % make model prediction binary:
   model = double(model>0);
   cc = corrcoef(model,rating);
   Corr_int(S) = cc(1,2)   
end

% display these results
figure(1); subplot(1,2,1)
plot(Corr_int,'sr');
hold on;
plot(Corr_saliency,'ok');
axis([0 11 0 0.8]);
bar([8,9],[mean(Corr_saliency),mean(Corr_int)]);
title('Figure 2 in paper');


% compute the average saliency differences dependent on subjects response
% To this end, concatenate all data to make a boxplot

alldataS = TOTAL{S}(:,4)';
alldataM = TOTAL{S}(:,5)';
for s=2:7
  alldataS = [alldataS TOTAL{S}(:,4)'];
  alldataM = [alldataM TOTAL{S}(:,5)'];
end
Ind = zeros(length(alldataM),1);
Ind(find(alldataS==0)) = 0; % 'equal' saliency
Ind(find(alldataS==2)) = 2; % sample 2 more salient
Ind(find(alldataS==1)) = 1; % sample 1

han  = subplot(1,2,2);
boxplot(alldataM,Ind,'notch','on','Color','kkk','labels',{'equal','stim1','stim2'});
set(gca,'YTick',[-0.8:0.4:0.8]);
axis([0.5 3.5 -0.8 0.8])
xlabel('Subjects choice');
ylabel('Model difference');
grid on



This package contains code, example sounds and data from the auditory saliency project,
published in  
"Mechanisms for allocating auditory attention: An auditory saliency map." 
Kayser C, Petkov C, Lippert M, Logothetis N. Current Biology 15(21), 1943-1947, 2005.


The code is made available to enhance developments around the ideas of auditory saliency. 
Direct use of this code should be acknowledged by citing the above paper. 



Saliency map functions:
-----------------------

Saliency_map.m: implements the saliency map (with sub-functions in the same file)

Example_saliency.m: illustrates how to use this using a simple example


Example behavioral  data:
--------------------------

Tutorial_saliencyexp1.m: reproduces a basic analysis of behavioral data from the paper. 

This exploits the actual data included in Data_saliencyexp_exp1.mat.


Example sounds in the study:
----------------------------

Sounds_saliencyexp_1.mat: provides the sound stimuli used in the above experiment. 



% load('Sounds_saliencyexp_1.mat','Sound')
% Sound is a 52x3 cell array with the sound samples, consisting of 52 sounds
% presented on 3 different backgrounds, i.e. Sound{sample,background}.
% During the experiments pairs of randomly chosen sounds were presented and subjects
% had to indicate which one was more salient. In between catch trials were presented,
% where the same sound was presented twice.

function out = frequency(in)
% This function obtains the auditory spectrum, simulating the information
% processing process of the early auditory system.
%
% INPUT
% -- in: Audio signal vector.Audio signal vector.

paras = [20,128,-2,0];
y = wav2aud(in, paras,'p_o');
ksize=10;
ksigma=fspecial('gaussian', ksize, 6);
out=conv2(y, ksigma, 'same');
out=out';
return



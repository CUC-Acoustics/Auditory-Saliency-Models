function [cr_r,cr_s] = rs_feature(in)
% This function is used to compute the rate and scale features.
%
% INPUT
% -- in: Audio signal vector.Audio signal vector.
paras = [20,128,-2,0];
y = wav2aud(in, paras,'p_o');
rv=2.^(0:0.5:8);
sv=2.^(-2:0.4:4);
cr = aud2cor(y, paras, rv, sv, 'kaya.cor');
cr_rs = mean(abs(cr), 4); 
cr_rs=squeeze(cr_rs);
cr_r=mean(cr_rs,1);
cr_r=squeeze(cr_r);
cr_s=mean(cr_rs,2);
cr_s=squeeze(cr_s);
return
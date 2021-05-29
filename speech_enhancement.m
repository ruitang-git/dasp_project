%% single channel
[bn,fs] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\babble_noise.wav");
[ssn,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\Speech_shaped_noise.wav");
[cs2,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech_2.wav");
[cs,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech.wav");
[ann,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\aritificial_nonstat_noise.wav");
seconds=10; range=fs*10+1:fs*(10+seconds);
bn=bn(range); ssn=ssn(range); cs2=cs2(range);cs=cs(range); ann=ann(range);
cs_bn = cs+bn;
cs_ssn = cs+ssn;
cs_ann = cs+ann;
cs2_bn = cs2+bn;
cs2_ssn = cs2+ssn;
cs2_ann = cs2+ann;
%% gain function
%% noise psd estimation
% VAD
% MMSE
time_int = 0.02;% second
x = cs_ssn;
sigma_n2 = mmse_based_noise_psd(x,fs,time_int);
%% signal psd estimation
%% multiple channel
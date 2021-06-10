clc;
clear;
%% single channel
[bn,fs] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\babble_noise.wav");
[ssn,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\Speech_shaped_noise.wav");
[cs2,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech_2.wav");
[cs,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech.wav");
[ann,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\aritificial_nonstat_noise.wav");
seconds=20; range=fs*3+1:fs*(3+seconds);
bn=bn(range); ssn=ssn(range); cs2=cs2(range);cs=cs(range); ann=ann(range);
%% assign data
signal = cs;
% noise = 0.2*ssn; % assgin noise
% noisy_signal = signal+noise;
% gaussian
noisy_signal = awgn(signal,35);
noise = noisy_signal-signal;
%% pure noise
% x = bn; % assgin noise
% time_int = 0.02;% second
% Nn = length(x);
% N_seg = time_int * fs;
% seg = ceil(Nn/N_seg); % number of segamentations
% zeropad = seg*N_seg -Nn;
% x_pad = [x',zeros(1,zeropad)];
% x_seg = reshape(x_pad,[N_seg,seg]); % one column = one segamentation
% pnn = periodogram(x_seg,rectwin(length(x_seg(:,1))),length(x_seg(:,1)));
% plot(pnn);
%pnn = periodogram(x_seg,rectwin(length(x_seg(:,1))),length(x_seg(:,1)),fs);
%% pure noise(overlap)
%% peridogam
L = 320; R = 160;
zero_pad = mod(length(noise),R);
x_pad = [noise;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
pnn = 2*pi*L*periodogram(x_seg,rectwin(length(x_seg(:,1))),length(x_seg(:,1)),'twosided');
%% smoothed periodogram
M = 10; % window size
T_sm = (M*R+R)/fs; %second
alpha = (T_sm*fs/R-1)/(T_sm*fs/R+1);
[~,col] = size(pnn);
pnn_s = pnn;
for j = M:col
    window = pnn(:,j-M+1:j);
    pnn_s(:,j) = window(:,1);
    for k = 2:M
        pnn_s(:,j) = alpha*pnn_s(:,j)+(1-alpha)*window(:,k);
    end
end
%% noise psd estimation(noisy signal)
%% VAD
%% MS
sigma_n2_ms = ms_based_noise_psd(noisy_signal,fs);
% figure;plot(real((sum(pnn))));hold on;plot(real((sum(pnn_s))));hold on;plot(real((sum(sigma_n2_ms))));
%% MMSE
sigma_n2_mmse = mmse_based_noise_psd(noisy_signal);
k = 25; % observe #k frequency
figure;plot(10*log(real(pnn(k,:))),'--','Linewidth',0.3);hold on;plot(10*log(real(pnn_s(k,:))));
hold on;plot(10*log(real(sigma_n2_ms(k,:))),'Linewidth',1.5);hold on;plot(10*log(real(sigma_n2_mmse(k,:))),'Linewidth',1.5);
title(['Noise PSD Estimation(K=',num2str(k),')']);
xlabel("time period");
ylabel("noise psd(dB)");
legend("periodogram","smoothed peridogram","minimum statistics(MS)","MMSE-SPP");
%% signal psd estimation
% method = 'DD';
% sigma_n2 = sigma_n2_mmse;
% sigma_s2 = signal_psd(noisy_signal,sigma_n2,method);
%% Wiener filter
alpha = 0.98;
sigma_n2 = pnn_s;
signal_est_rec = wiener(noisy_signal,sigma_n2,alpha);
%% multiple channel
N = 100;
figure;
subplot(3,1,1);
plot(noisy_signal(N:end));
subplot(3,1,2);
plot(signal(N:end));
subplot(3,1,3);
plot(signal_est_rec(N:end));
ylim([-0.5 0.5])
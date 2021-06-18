clc;
clear;
%% single channel
[bn,fs] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\babble_noise.wav");
[ssn,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\Speech_shaped_noise.wav");
[cs2,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech_2.wav");
[cs,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\clean_speech.wav");
[ann,~] = audioread("D:\slides\EE4182 Digital Audio and Speech Processing\project\sound files for mini-project\aritificial_nonstat_noise.wav");
seconds=10; range=fs*3+1:fs*(3+seconds);
bn=bn(range); ssn=ssn(range); cs2=cs2(range);cs=cs(range); ann=ann(range);
%% assign data
signal = cs;
noise = 0.2*bn; % assgin noise
noisy_signal = signal+noise;
% gaussian
% noisy_signal = awgn(signal,35);
% noise = noisy_signal-signal;
%% assign parameters
K = 10;
alpha = 0.98;
method = 'DD';
% method = 'ML';
% method = 'DD';
noise_psd_estimator = 'MMSE';
% noise_psd_estimator = 'MS';
% noise_psd_estimator = 'MMSE';
gain_function = 'spectral_subtrac';
% gain_function = 'wiener';
% gain_function = 'spectral_subtrac';


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
alpha_sp = (T_sm*fs/R-1)/(T_sm*fs/R+1);
[~,col] = size(pnn);
pnn_s = pnn;
for j = M:col
    window = pnn(:,j-M+1:j);
    pnn_s(:,j) = window(:,1);
    for k = 2:M
        pnn_s(:,j) = alpha_sp*pnn_s(:,j)+(1-alpha_sp)*window(:,k);
    end
end
%% pure signal(overlap)
% %% peridogam
% zero_pad = mod(length(signal),R);
% x_pad = [signal;zeros(zero_pad,1)];
% x_seg = [];
% i = 1;
% while i~=length(x_pad)/R
%     x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
%     i = i+1;
% end
% pss = 2*pi*L*periodogram(x_seg,rectwin(length(x_seg(:,1))),length(x_seg(:,1)),'twosided');
% %% smoothed periodogram
% pss_s = pss;
% for j = M:col
%     window = pss(:,j-M+1:j);
%     pss_s(:,j) = window(:,1);
%     for k = 2:M
%         pss_s(:,j) = alpha_sp*pss_s(:,j)+(1-alpha_sp)*window(:,k);
%     end
% end
%% noise psd estimation(noisy signal)
%% VAD
%% MS
sigma_n2_ms = ms_based_noise_psd(noisy_signal,fs);
% figure;plot(real((sum(pnn))));hold on;plot(real((sum(pnn_s))));hold on;plot(real((sum(sigma_n2_ms))));
%% MMSE
sigma_n2_mmse = mmse_based_noise_psd(noisy_signal);
%% signal psd estimation
% method = 'DD';
% sigma_n2 = sigma_n2_mmse;
% sigma_s2 = signal_psd(noisy_signal,sigma_n2,method);
%% Wiener filter
if strcmp(noise_psd_estimator,"MMSE")
    sigma_n2 = sigma_n2_mmse;
elseif strcmp(noise_psd_estimator,"MS")
    sigma_n2 = sigma_n2_ms;
end   
if strcmp(gain_function,'wiener')
    [signal_est_rec,sigma_s2] = wiener(noisy_signal,sigma_n2,alpha,K,method);
elseif strcmp(gain_function,'spectral_subtrac')
    [signal_est_rec,sigma_s2] = spectral_subtraction(noisy_signal,sigma_n2,K);
end
%% performance assessment
evaluate_denoising_metrics(signal,signal_est_rec);
%% plot noise psd
k = 25; % observe #k frequency
figure;plot(10*log(real(pnn(k,:))),'--','Linewidth',0.3);...
    hold on;plot(10*log(real(pnn_s(k,:))));...
    hold on;plot(10*log(real(sigma_n2_ms(k,:))),'Linewidth',1.5);...
    hold on;plot(10*log(real(sigma_n2_mmse(k,:))),'Linewidth',1.5);
title(['Noise PSD Estimation(K=',num2str(k),')']);
xlabel("time period");
ylabel("noise psd(dB)");
legend("periodogram","smoothed peridogram","minimum statistics(MS)","MMSE-SPP");
%% plot signal psd
% figure;plot(10*log(real(pss(k,:))),'--','Linewidth',0.5);...
%     hold on;plot(10*log(real(pnn_s(k,:))));...
%     %hold on;plot(10*log(real(sigma_s2(k,:))),'Linewidth',1);
%% plot denoised signal
N = 5e3;
figure;
subplot(3,1,1);
plot(noisy_signal(N:end));
subplot(3,1,2);
plot(signal(N:end));
subplot(3,1,3);
plot(signal_est_rec(N:end));
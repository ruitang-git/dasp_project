%% signal psd estimation
function sigma_s2 = signal_psd(x,sigma_n2,method)
% peridogam
L = 320; R = 160;
zero_pad = mod(length(x),R);
x_pad = [x;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
pss = periodogram(x_seg,rectwin(length(x_seg(:,1))),length(x_seg(:,1)));
if strcmp(method,'bartlett')
% Bartlett
K = 10; % window size
[~,col] = size(pss);
pss_s = pss;
for j = K:col
    window = pss(:,j-K+1:j);
    pss_s(:,j) = sum(window,2)/K;
end
sigma_s2 = pss_s - sigma_n2;
elseif strcmp(method,'DD')
alpha = 0.97;
sigma_s2 = pss - sigma_n2;
for i = 2:size(sigma_s2,2)
    sigma_s2(:,i) = alpha*(sigma_s2(:,i-1)) + (1-alpha)*max(sigma_s2(:,i),0);
end
end
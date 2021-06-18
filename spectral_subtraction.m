function [signal_est_rec,sigma_s2] = spectral_subtraction(x,sigma_n2,K)
L = 320; R = 160;
zero_pad = mod(length(x),R);
x_pad = [x;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
% window
% win_seg = hann(L);
% window =[];
% for i = 1:size(x_seg,2)
%     window(:,i) = win_seg;
% end
% x_seg = x_seg.*window;

y = fft(x_seg);
signal_est_rec = [];
sigma_s2 = [];
% y = y(1:R+1,:);

for i = 1:size(y,2)
    if i>=K
        %snr = max(sum((abs(y(:,i-K+1:i)).^2),2)./(K*sigma_n2(:,i))-1,0);
        y_sqr = sum((abs(y(:,i-K+1:i)).^2),2)./K;
    else
        %snr = max((abs(y(:,i)).^2)./sigma_n2(:,i)-1,0);
        y_sqr = abs(y(:,i)).^2;
    end
    buffer = max(y_sqr - sigma_n2(:,i),0);
    sigma_s2 = [sigma_s2,buffer];
    s_est = sqrt(buffer);
    % s_est_rec = [s_est;flip(s_est(2:end-1,:))];
    signal_est = ifft(s_est);
    if i == 1
        signal_est = real(signal_est);
    else
        signal_est = real(signal_est(L-R+1:end));
    end
    signal_est_rec = [signal_est_rec;signal_est];
end

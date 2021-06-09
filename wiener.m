function signal_est_rec = wiener(x,sigma_n2,sigma_s2)
H = sigma_s2./(sigma_n2+sigma_s2);
L = 320; R = 160;
zero_pad = mod(length(x),R);
x_pad = [x;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
% fft
y = fft(x_seg);
y = y(1:R+1,:);
s_est = real(H).*y;
s_est_rec = [s_est;flip(s_est(2:end-1,:))];
signal_est = ifft(s_est_rec);
signal_est_rec = signal_est(:,1);
for i = 2:size(signal_est,2)
    seg = signal_est(L-R+1:end,i);
    signal_est_rec = [signal_est_rec;seg];
end
signal_est_rec = real(signal_est_rec);

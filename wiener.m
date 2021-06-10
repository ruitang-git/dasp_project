function signal_est_rec = wiener(x,sigma_n2,alpha)
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
% win_seg = hann(L);
% window =[];
% for i = 1:size(x_seg,2)
%     window(:,i) = win_seg;
% end
% x_seg = x_seg.*window;

y = fft(x_seg);
% y = y(1:R+1,:);
s_est = zeros(size(y,1),1);
signal_est_rec = [];
for i = 1:size(y,2)

    snr = alpha * abs(s_est).^2 ./ sigma_n2(:,i) +...
        (1-alpha) * max((abs(y(:,i)).^2)./sigma_n2(:,i)-1,0);

    s_est = (snr./(1+snr)).* y(:,i);
    % s_est_rec = [s_est;flip(s_est(2:end-1,:))];
    signal_est = ifft(s_est);
    if i == 1
        signal_est = signal_est;
    else
        signal_est = signal_est(L-R+1:end);
    end
    signal_est_rec = [signal_est_rec;signal_est];
end
signal_est_rec = real(signal_est_rec);

function [signal_est_rec,sigma_s2] = wiener(x,sigma_n2,alpha,K,method)
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
if strcmp(method,"DD")
s_est = zeros(size(y,1),1);
    for i = 1:size(y,2)
        snr = alpha * abs(s_est).^2 ./ sigma_n2(:,i) +...
            (1-alpha) * max((abs(y(:,i)).^2)./sigma_n2(:,i)-1,0);
        buffer = snr.*sigma_n2(:,i);
        sigma_s2 = [sigma_s2,buffer];
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
elseif strcmp(method,"ML")
    for i = 1:size(y,2)
        if i>=K
            %snr = max(sum((abs(y(:,i-K+1:i)).^2),2)./(K*sigma_n2(:,i))-1,0);
            snr = sum((abs(y(:,i-K+1:i)).^2),2)./(K*sigma_n2(:,i))-1;
        else
            %snr = max((abs(y(:,i)).^2)./sigma_n2(:,i)-1,0);
            snr = abs(y(:,i)).^2./sigma_n2(:,i)-1;
        end
        buffer = snr.*sigma_n2(:,i);
        sigma_s2 = [sigma_s2,buffer];
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
%     y2 = abs(y).^2;
%     y2_s = y2;
%     for j = K:size(y2,2)
%         window = y2(:,j-K+1:j);
%         y2_s(:,j) = sum(window,2)/K;
%     end
%     sigma_s2 = max(y2_s - sigma_n2,0);
%     H = sigma_s2./(sigma_s2+sigma_n2);
%     signal_est = ifft((sigma_s2./(sigma_s2+sigma_n2)).*y);
%     for i = 1:size(y2,2)
%         if i == 1
%             buffer = signal_est(:,i);
%         else
%             buffer = signal_est(L-R+1:end,i);
%         end
%         signal_est_rec = [signal_est_rec;buffer];
%     end
end
    


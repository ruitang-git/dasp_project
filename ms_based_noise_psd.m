%% MS based noise psd estimation
function sigma_n2 = ms_based_noise_psd(x,fs)
L = 320; R = 160;
zero_pad = mod(length(x),R);
x_pad = [x;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
M = 10;
T_sm = (M*R+R)/fs; %second
alpha = (T_sm*fs/R-1)/(T_sm*fs/R+1);
f = fft(x_seg);
f = f(1:R+1,:);
f_2 = abs(f).^2;
p = f_2;
[row,col] = size(p);
p_s = p;
for j = M:col
    window = p(:,j-M+1:j);
    p_s(:,j) = window(:,1);
    for k = 2:M
        p_s(:,j) = alpha*p_s(:,j)+(1-alpha)*window(:,k);
    end
end
sigma_n2 = p_s;
N = 50;
for k = 1:N-1
    for l = 1:row
        sigma_n2(l,k) = min(p_s(l,1:k));
    end
end
for k = N:col
    for l = 1:row
        sigma_n2(l,k) = min(p_s(l,k-N+1:k));
    end
end
sigma_n2 = sigma_n2./(2*pi*L);
sigma_n2(2:end-1,:) = 2*sigma_n2(2:end-1,:);



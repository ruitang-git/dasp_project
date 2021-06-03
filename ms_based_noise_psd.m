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
M = 20;
T_sm = (M*R+R)/fs; %second
alpha = (T_sm*fs/R-1)/(T_sm*fs/R+1);
f = fft(x_seg);
f = f(1:R+1,:);
f_2 = abs(f).^2;
p = f_2;
[row,col] = size(p);
for j = 2:col
    p(:,j) = alpha.*p(:,j-1)+(1-alpha)*p(:,j);
end
sigma_n2 = p;
for k = M:col
    for l = 1:row
        sigma_n2(l,k) = min(p(l,k-M+1:k));
    end
end
sigma_n2 = sigma_n2./(2*pi*L);
sigma_n2(2:end-1,:) = 2*sigma_n2(2:end-1,:);



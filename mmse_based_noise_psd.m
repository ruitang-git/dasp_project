%% MMSE
function sigma_n2 = mmse_based_noise_psd(x)
% segamentation
L = 320; R = 160;
zero_pad = mod(length(x),R);
x_pad = [x;zeros(zero_pad,1)];
x_seg = [];
i = 1;
while i~=length(x_pad)/R
    x_seg(:,i) = x_pad(i*R-R+1:(i+1)*R);
    i = i+1;
end
[row,col]= size(x_seg);
% fft
f = fft(x_seg);
% f = f(1:R+1,:);
% psd estimation
p_H0 = 0.5;p_H1 = 0.5;
alpha = 0.8;
noise_psd = [];
sigma_n2 = [];
cosin = 10^1.5;
for j = 1:5
    %for k = 1:R+1
    for k = 1:L
       noise_psd(k,j) = f(k,j)*f(k,j)';
       if j == 1
           sigma_n2(k,j) = noise_psd(k,j);
       else
           sigma_n2(k,j) = alpha*sigma_n2(k,j-1)+(1-alpha)*noise_psd(k,j);   
       end
    end
end
for j = 6:col
    % for k = 1:R+1
    for k = 1:L
        p_H1y = (1+(p_H0/p_H1)*(1+cosin)*exp(-(f(k,j)*f(k,j)'*cosin)/((1+cosin)*sigma_n2(k,j-1))))\1;
        p_H0y = 1-p_H1y;
        noise_psd(k,j) = p_H0y*f(k,j)*f(k,j)'+ p_H1y*sigma_n2(k,j-1);
        sigma_n2(k,j) = alpha*sigma_n2(k,j-1)+(1-alpha)*noise_psd(k,j);
    end
end
% sigma_n2 = sigma_n2./(2*pi*L);
% sigma_n2(2:end-1,:) = 2*sigma_n2(2:end-1,:);


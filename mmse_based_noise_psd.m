%% MMSE
function sigma_n2 = mmse_based_noise_psd(x,fs,time_int)
% segamentation
N = length(x);
N_seg = time_int * fs;
seg = ceil(N/N_seg); % number of segamentations
zeropad = seg*N_seg -N;
x_pad = [x',zeros(1,zeropad)];
x_seg = reshape(x_pad,[N_seg,seg]); % one column = one segamentation
% fft
f = [];
for j = 1:seg
    f(:,j) = fft(x_seg(:,1));
end 
% psd estimation
p_H0 = 0.5;p_H1 = 0.5;
alpha = 0.8;
noise_psd = [];
sigma_n2 = [];
cosin = 10^1.5;
for j = 1:5
    for k = 1:N_seg
       noise_psd(k,j) = f(k,j)*f(k,j)';
       if j == 1
           sigma_n2(k,j) = noise_psd(k,j);
       else
           sigma_n2(k,j) = alpha*sigma_n2(k,j-1)+(1-alpha)*noise_psd(k,j);   
       end
    end
end
for j = 6:seg
    for k = 1:N_seg
        p_H1y = (1+p_H0/p_H1*(1+cosin)*exp(-(f(k,j)*f(k,j)'*cosin)/((1+cosin)*sigma_n2(k,j-1))));
        p_H0y = 1-p_H1y;
        noise_psd(k,j) = p_H0y*f(k,j)*f(k,j)'+ p_H1y*sigma_n2(k,j-1);
        sigma_n2(k,j) = alpha*sigma_n2(k,j-1)+(1-alpha)*noise_psd(k,j);
    end
end
sigma_n2 = sigma_n2./N_seg;

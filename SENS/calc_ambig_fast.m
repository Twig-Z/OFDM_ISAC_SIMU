function [chi_zero_delay, chi_zero_doppler] = calc_ambig_fast(s_t, tau_vec, fd_vec, fs)

s_t = reshape(s_t, [], 1);
N_total = length(s_t);
t = (0:N_total-1).' / fs;

N_fd = length(fd_vec);

%% 零多普勒切片 χ(τ,0) - 直接用相关函数
% χ(τ,0) = |∫s(t)s*(t-τ)dt|^2 = 自相关函数
[acf, lags] = xcorr(s_t, s_t);
acf_power = abs(acf);  % 去掉平方符合你的习惯

% 插值到tau_vec
lags_time = lags / fs;
chi_zero_doppler = interp1(lags_time, acf_power, tau_vec, 'spline', 0);
chi_zero_doppler = chi_zero_doppler(:);

%% 零时延切片 χ(0,fd) - 直接用FFT
% χ(0,fd) = |∫|s(t)|^2 * exp(-j2πfd*t)dt|^2
% 即|s(t)|^2的傅里叶变换
s_power = abs(s_t).^2;  % [N_total x 1]
% fd_vec [1 x N_fd], t [N_total x 1]
exp_matrix = exp(-1j * 2 * pi * t * fd_vec(:).');  % [N_total x N_fd]
chi_zero_delay = abs(s_power.' * exp_matrix) / fs;  % [1 x N_fd]
chi_zero_delay = chi_zero_delay(:);
    
% 只归一化切片
max_val = max(max(chi_zero_doppler), max(chi_zero_delay));  
chi_zero_delay = chi_zero_delay / max_val;  
chi_zero_doppler = chi_zero_doppler / max_val;
end

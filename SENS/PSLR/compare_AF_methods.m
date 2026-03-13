

clc; clear; close all;
%% ========================================================================
%  模糊函数计算方法对比实验
%
%  功能：
%  本程序用于对比三种不同的模糊函数（Ambiguity Function, AF）计算方法，
%  并分析其计算精度与计算复杂度之间的差异。程序生成OFDM信号后，
%  分别使用三种方法计算模糊函数或其切片，并对结果进行可视化比较。
%
% ------------------------------------------------------------------------
%  方法1：calc_ambig （标准模糊函数计算）
%
%  - 完全按照模糊函数的连续定义公式进行离散计算：
%
%      χ(τ, f_d) = | ∫ s(t) · s*(t-τ) · e^{-j2πf_d t} dt |
%
%  - 对每一个 (τ, f_d) 点进行逐点积分计算
%  - 计算量非常大，复杂度约为：
%
%      O(N_signal × N_tau × N_fd)
%
%  - 优点：
%        • 计算结果最精确
%        • 与理论定义完全一致
%  - 缺点：
%        • 计算速度非常慢
%
%  因此该方法作为“基准方法（Ground Truth）”用于结果对比。
%
% ------------------------------------------------------------------------
%  方法2：calc_ambig_xcorr2 （基于 xcorr2 的近似计算）
%
%  - 利用二维互相关 xcorr2 对模糊函数进行近似计算
%  - 本质上等价于：
%
%      χ(τ,0) ≈ xcorr(s(t), s(t))
%
%  - 对零多普勒切片 χ(τ,0) 的结果与方法1基本一致
%  - 但对于零时延切片 χ(0,f_d) 的结果存在一定偏差
%
%  - 优点：
%        • 计算速度显著提高
%        • 实现简单
%
%  - 缺点：
%        • 本质是近似方法
%        • 对完整模糊函数平面精度有限
%
% ------------------------------------------------------------------------
%  方法3：calc_ambig_slice （快速切片计算方法）
%
%  为后续快速模糊函数计算设计，仅计算两个关键切片：
%
%      1）零多普勒切片 χ(τ,0)
%      2）零时延切片 χ(0,f_d)
%
%  具体实现思路：
%
%  • 零多普勒切片 χ(τ,0)
%    通过向量化分数时延（fractional delay）实现：
%
%        s(t-τ) 使用 interp1 进行分数采样延迟
%
%    并一次性构造所有 τ 的查询矩阵，实现完全向量化计算。
%
%  • 零时延切片 χ(0,f_d)
%
%        χ(0,f_d) = | ∫ |s(t)|² e^{-j2πf_d t} dt |
%
%    本质上等价于对信号功率 |s(t)|² 进行频谱分析，
%    通过构造指数矩阵实现向量化计算。
%
%  - 优点：
%        • 计算速度极快
%        • 适合实时或大规模计算
%  - 缺点：
%        • 仅计算两个切片，不得到完整 AF 平面
%
% ------------------------------------------------------------------------
%  实验输出：
%
%  程序将输出以下对比结果：
%
%  1. 完整模糊函数 χ(τ,f_d) 平面对比（方法1 vs 方法2）
%  2. 模糊函数三维图
%  3. 两种方法的差值图
%  4. 三种方法的切片对比：
%        • Zero-Doppler 切片 χ(τ,0)
%        • Zero-Delay 切片 χ(0,f_d)
%
%  用于评估不同算法的精度与计算效率。
%
% ------------------------------------------------------------------------
%  作者：ZZZ
%  用途：OFDM感知信号模糊函数计算方法研究
% ========================================================================
%% 规模
flag = 1;
switch flag
    case 1
        n_fft = 512;
        n_sc  = 200;
        n_sym = 40;
    case 2
        n_fft = 1024;
        n_sc  = 600;
        n_sym = 100;
    case 3
        n_fft = 2048;
        n_sc  = 600;
        n_sym = 100;
end

%% 信号参数
fc      = 10e9;
delta_f = 120e3;

T_O  = 1 / delta_f;
T_cp = 0.58e-6;
T_sym = T_O + T_cp;

fs = n_fft * delta_f;

n_cp = round(T_cp * fs);

%% 延迟多普勒范围
c = 3e8;
lambda = c / fc;

max_range = 100;
tau_max = max_range / c;

N_tau = 201;
tau_vec = linspace(-tau_max, tau_max, N_tau);

max_velocity = 100;
fd_max = 2 * max_velocity / lambda;

N_fd = 201;
fd_vec = linspace(-fd_max, fd_max, N_fd);

%% 随机16QAM信号
qam_symbols = qammod(randi([0 15], n_sc, n_sym), 16, 'UnitAveragePower', true);

S_fd1 = zeros(n_fft, n_sym);

sc_start = 2;
sc_end   = n_sc + 1;

S_fd1(sc_start:sc_end, :) = qam_symbols;

s_t = reshape(ifft(S_fd1, n_fft), [], 1);

%% 计算模糊函数
fprintf('\n计算模糊函数...\n');

tic
[chi1, chi_0delay1, chi_0doppler1] = calc_ambig(s_t, tau_vec, fd_vec, fs);
toc

tic
[chi2, chi_0delay2, chi_0doppler2] = calc_ambig_xcorr2(s_t, tau_vec, fd_vec, fs);
toc

tic
[chi_0delay3, chi_0doppler3] = calc_ambig_slice(s_t, tau_vec, fd_vec, fs);
toc

%% 完整模糊函数
figure(1);

subplot(1,2,1)
imagesc(fd_vec, tau_vec, chi1)
set(gca,'YDir','normal')
xlabel('多普勒 f_d (Hz)')
ylabel('时延 τ (s)')
title('模糊函数 χ(τ, f_d) 平面trapz积分法')
colorbar

subplot(1,2,2)
imagesc(fd_vec, tau_vec, chi2)
set(gca,'YDir','normal')
xlabel('多普勒 f_d (Hz)')
ylabel('时延 τ (s)')
title('模糊函数 χ(τ, f_d) 平面xcorr2法')
colorbar

figure(2)

subplot(1,2,1)
mesh(tau_vec, fd_vec, chi1')
xlabel('时延 τ (s)')
ylabel('多普勒 f_d (Hz)')
zlabel('幅度')
title('模糊函数χ(τ, f_d) 立体trapz积分')
axis tight
view(45,30)

subplot(1,2,2)
mesh(tau_vec, fd_vec, chi2')
xlabel('时延 τ (s)')
ylabel('多普勒 f_d (Hz)')
zlabel('幅度')
title('模糊函数χ(τ, f_d) 立体xcorr2法')
axis tight
view(45,30)

%% 差值
figure(3)

subplot(1,2,1)
imagesc(fd_vec, tau_vec, abs(chi1-chi2))
set(gca,'YDir','normal')
xlabel('多普勒 f_d (Hz)')
ylabel('时延 τ (s)')
title('模糊函数 χ(τ, f_d) 平面差值')
colorbar

subplot(1,2,2)
mesh(tau_vec, fd_vec, abs(chi1'-chi2'))
xlabel('时延 τ (s)')
ylabel('多普勒 f_d (Hz)')
zlabel('幅度')
title('模糊函数χ(τ, f_d) 立体差值')
axis tight
view(45,30)

%% 切片
figure(4)

subplot(2,2,1)
plot(tau_vec, chi_0doppler1 + 1e-12,'g','LineWidth',1.5)
hold on
plot(tau_vec, chi_0doppler2 + 1e-12,'r','LineWidth',1.5)
plot(tau_vec, chi_0doppler3 + 1e-12,'b','LineWidth',1.5)

xlabel('时延 τ (s)')
ylabel('幅度')
title('Zero-Doppler切片 χ(τ,0)对比')
legend('信号1','信号2','信号3','Location','best')
grid on
hold off

subplot(2,2,2)
plot(tau_vec, abs(chi_0doppler1-chi_0doppler2),'g','LineWidth',1.5)
hold on
plot(tau_vec, abs(chi_0doppler1-chi_0doppler3),'r','LineWidth',1.5)

xlabel('时延 τ (s)')
ylabel('幅度')
title('Zero-Doppler切片 χ(τ,0)差值')
legend('信号1-信号2','信号1-信号3','Location','best')
grid on
hold off

subplot(2,2,3)
plot(fd_vec, chi_0delay1 + 1e-12,'g','LineWidth',1.5)
hold on
plot(fd_vec, chi_0delay2 + 1e-12,'r','LineWidth',1.5)
plot(fd_vec, chi_0delay3 + 1e-12,'b','LineWidth',1.5)

xlabel('多普勒 f_d (Hz)')
ylabel('幅度')
title('Zero-Delay切片 χ(0,f_d)对比')
legend('信号1','信号2','信号3','Location','best')
grid on

subplot(2,2,4)
plot(fd_vec, abs(chi_0delay1-chi_0delay2) + 1e-12,'g','LineWidth',1.5)
hold on

chi222 = chi_0delay3';

plot(fd_vec, abs(chi_0delay1-chi_0delay3) + 1e-12,'r','LineWidth',1.5)

xlabel('多普勒 f_d (Hz)')
ylabel('幅度')
title('Zero-Delay切片 χ(0,f_d)差值')
legend('信号1-信号2','信号1-信号3','Location','best')
grid on
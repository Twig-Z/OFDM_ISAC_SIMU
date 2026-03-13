%%
clc;
clear;
close all;
tic;
disp("simulation start");

SystemSetup_simu
currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_ttt.mat'));


SNR=0:2:40;                     % 仿真信噪比
sim_num = 5;                    % 仿真次数（用于算平均误码率）

%% 瑞利信道设计
maxDoppler = 300;                                          %多普勒频偏
pathPower = [-1.0 -1.0 -1.0 0 0 0 -3.0 -5.0 -7.0];
pathDelays = [0 50 120 200 230 500 1600 2300 5000]*1e-9;

rchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
    'MaximumDopplerShift',maxDoppler,'FadingTechnique','Sum of sinusoids');

rng(1);

tx_data = randi([0 1], 1, N_b_data);
data_bit = tx_data;
[txSig,~] = OFDM_proc_tx(tx_data);

%% 信道（通过多经瑞利信道、AWGN信道）
Ber1 = zeros(1,length(SNR)); % 译码前误码率统计
Ber2 = zeros(1,length(SNR)); % 译码后误码率统计

for snr_ind=1:length(SNR) 
    fprintf('正在仿真SNR = %d dB...\n', SNR(snr_ind));
    
    for sim_ind = 1:sim_num    
        % 信道传输
        channel_out = step(rchan, txSig);
        rxSig = awgn(channel_out, SNR(snr_ind), 'measured');

[rx_data,y_stream,~] = OFDM_proc_rx(rxSig);
%% 计算误码率
        %[~, ber1] = biterr(y_de_bit(1:length(data_intlv)), data_intlv);
        [~, ber1] = biterr(rx_data(1:length(data_bit)), data_bit);
        Ber1(snr_ind) = Ber1(snr_ind) + ber1;
        [~, ber2] = biterr(rx_data(1:length(data_bit)), data_bit);
        Ber2(snr_ind) = Ber2(snr_ind) + ber2;
    end
    
    Ber1(snr_ind) = Ber1(snr_ind) / sim_num;
    Ber2(snr_ind) = Ber2(snr_ind) / sim_num;
    if Ber2(snr_ind) == 0
        Ber1(snr_ind) = NaN;   % 用NaN代替，图上自动断开
        Ber2(snr_ind) = NaN;
    end  
    fprintf("SNR = %d dB, 译码前BER = %.5f, 译码后BER = %.5f\n", SNR(snr_ind), Ber1(snr_ind), Ber2(snr_ind));
end


figure(2);
semilogy(SNR, Ber1, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
semilogy(SNR, Ber2, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
SNR_lin = 10.^(SNR/10);
BER_theory = qfunc(sqrt(2*SNR_lin));   % AWGN下QPSK理论值
BER_theory(BER_theory < 1e-5) = NaN;  % 低于1e-5不画
semilogy(SNR, BER_theory, 'k--', 'LineWidth', 1.5);
legend('解交织、卷积码译码前', '解交织、卷积码译码后', 'AWGN理论最小值','Location', 'best');
xlabel('SNR (dB)');
ylabel('BER');
title('瑞利衰落信道下误比特率曲线');
xlim([min(SNR), max(SNR)]);


figure(3);
subplot(2,1,1);
stem(0:30, data_bit(1:31));
ylabel('Amplitude');
title('发送数据（前31个比特）');
grid on;

subplot(2,1,2);
stem(0:30, rx_data(1:31));
ylabel('Amplitude');
title('接收数据（前31个比特）');
grid on;

% 显示接收信号星座图

scatterplot(y_stream);
title('接收信号星座图（均衡后）');

fprintf('仿真完成，总耗时: %.2f 秒\n', toc);
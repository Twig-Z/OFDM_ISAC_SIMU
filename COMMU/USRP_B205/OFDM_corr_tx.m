%% ==============================================================
%  OFDM_corr_tx.m
%  功能：
%  生成带扩频和信道编码的 OFDM 发射信号，并通过 USRP 发送。
%
%  主要流程：
%  1. 根据系统带宽自动计算 OFDM 参数
%  2. 生成随机比特
%  3. 卷积编码
%  4. QPSK 调制
%  5. m序列扩频
%  6. 插入导频
%  7. OFDM 子载波映射
%  8. IFFT + 循环前缀
%  9. 添加 PSS 同步序列
% 10. USRP 连续发射
% ==============================================================

clc;
clear;
 
SystemSetup
currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_ttt.mat'));
%load('parameter_ttt.mat')
%% ================= 生成比特 =================

rng(1);

P_data = randi([0 1], 1, N_b_data);

%% ================= 卷积编码 =================

trellis = poly2trellis(L, [133 171]);

code_data = convenc(P_data, trellis);

%% ================= QPSK调制 =================

data_temp1 = reshape(code_data, log2(M), [])';

data_temp2 = bi2de(data_temp1);

modu_data = pskmod(data_temp2, M, pi/M);

%% ================= 扩频 =================

code = spread_mseq(stage, ptap1, regi1, data_row);

code = code * 2 - 1;

modu_data = reshape(modu_data, data_row, []);

spread_data = spread(modu_data, code);

spread_data = reshape(spread_data, [], 1);

modu_data = spread_data;

%% ================= 插入导频 =================

rng(1);

pilot_seq = 2 * (randi([0 1], pilot_num, data_col * spread_ratio)) - 1;

data = zeros(N_sc, data_col * spread_ratio);

data(P_f_station, :) = pilot_seq;

if data_row * data_col > length(modu_data)

    data2 = [modu_data; zeros(data_row * data_col - length(modu_data), 1)];

else

    data2 = modu_data;

end

data_seq = reshape(data2, data_row, data_col * spread_ratio);

data(data_station, :) = data_seq;

%% ================= 子载波映射 =================

data = [zeros(pos_row - 1, data_col * spread_ratio); data];

data3 = [data; zeros(N_fft - N_sc - pos_row + 1, data_col * spread_ratio)];

%% ================= OFDM调制 =================

ifft_data = ifft(data3);

Tx_cd = [ifft_data(N_fft - N_cp + 1:end, :); ifft_data];

%% ================= PSS同步 =================

PSS_code = spread_mseq(PSS_stage, PSS_ptap1, PSS_regi1, 1);

PSS_code2 = PSS_code * 2 - 1;

PSS = sqrt(1/N_fft) * PSS_code2(:);

%% ================= 串并转换 =================

Tx_data = reshape(Tx_cd, [], 1);

zeros_ahead = zeros(N_zeros, 1);

txSig = [zeros_ahead; PSS; Tx_data];

txSig = txSig / max(abs(txSig));


%% ================= USRP TX =================

radio = comm.SDRuTransmitter( ...
    'Platform', 'B200', ...
    'SerialNum', '2409004', ...
    'MasterClockRate', masterClock, ...
    'CenterFrequency', fc, ...
    'Gain', gain_tx, ...
    'InterpolationFactor', interp);

radio.ChannelMapping = 1;

radio.UnderrunOutputPort = true;

disp('TX running...');

%% ================= 发送 =================

while true

    radio(txSig);

    underrun = radio(txSig);

    if underrun
        disp('Underrun!');
    end

end

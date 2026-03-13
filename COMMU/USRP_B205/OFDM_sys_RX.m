%% ============================SNR_TIME_DOMAIN==================================
%  OFDM_sys_rx.m
%  功能：
%  接收带扩频与信道编码的 OFDM 信号，并完成解调与 BER 计算。
%
%  主要流程：
%  1. 根据系统带宽自动计算 OFDM 参数
%  2. USRP 接收信号
%  3. PSS 前导相关同步
%  4. OFDM 符号定位
%  5. 去循环前缀
%  6. FFT
%  7. 导频信道估计
%  8. 频域均衡
%  9. 解扩
% 10. QPSK 解调
% 11. 维特比译码
% 12. BER 与 SNR 统计
% ==============================================================

clc;
clear;
close all;

[fList, pList] = matlab.codetools.requiredFilesAndProducts('OFDM_sys_TX.m');

% 查看返回的自定义文件列表
disp('该文件调用的非系统函数/文件有：');
fList'
%% ================= System Setup =================

stage = 1;
ptap1 = [1];
regi1 = [1];

Nd = 14;

M = 4;

L = 7;
tblen = 6 * L;

P_f_inter = 3;

%% ================= OFDM系统参数 =================

B = 5e6;

mu = 0;
delta_f = 15625 * 2^mu;

eta = 0.8;

spread_ratio = (2^stage - 1);

B_leackage = 0.5e6;

B0 = (B - B_leackage) / spread_ratio;

N_sc = floor(eta * B0 / delta_f);

N_sc0 = floor((B0 + B_leackage) / delta_f);

N_fft = 2^ceil(log2(N_sc0));

N_cp = floor(0.07 * N_fft);

N_symbol = N_fft + N_cp;

pos_row = ceil(B_leackage / delta_f / 2);

%% ================= 子载波结构 =================

data_station = [];

for img = 1:N_sc
    if mod(img, P_f_inter) ~= 1
        data_station = [data_station img];
    end
end

data_row = length(data_station);

%% ================= 编码参数 =================

R_channelcoding = 1/2;

data_col = Nd;

P_f = 1/sqrt(2) + 1i/sqrt(2);

P_f_station = 1:P_f_inter:N_sc;

pilot_num = length(P_f_station);

rng(1);

pilot_seq = 2 * (randi([0 1], pilot_num, data_col * spread_ratio)) - 1;

%% ================= 扩频码 =================

code = spread_mseq(stage, ptap1, regi1, data_row);

code = code * 2 - 1;

%% ================= 数据产生 =================

N_RE_data = data_row * Nd;

N_b_data = N_RE_data * log2(M) * R_channelcoding;

rng(1);

P_data = randi([0 1], 1, N_b_data);

%% ================= 卷积编码 =================

trellis = poly2trellis(L, [133 171]);

code_data = convenc(P_data, trellis);

%% ================= QPSK调制 =================

data_temp1 = reshape(code_data, log2(M), [])';

data_temp2 = bi2de(data_temp1);

modu_data = pskmod(data_temp2, M, pi/M);

%% ================= PSS序列 =================

PSS_stage = 9;

PSS_ptap1 = [4 9];

PSS_regi1 = [1 1 0 1 1 0 1 1 1];

PSS_code = spread_mseq(PSS_stage, PSS_ptap1, PSS_regi1, 1);

PSS_code2 = PSS_code * 2 - 1;

PSS_transpose = PSS_code2.';

%% ================= USRP参数 =================

fc = 915e6;

gain_rx = 40;

samplesPerFrame = 64000;

masterClock = 32e6;

fs = delta_f * N_fft;

decim = round(masterClock / fs);

%% ================= USRP接收机 =================

radio = comm.SDRuReceiver( ...
    'Platform','B200', ...
    'SerialNum','2409028', ...
    'MasterClockRate',masterClock, ...
    'CenterFrequency',fc, ...
    'Gain',gain_rx, ...
    'DecimationFactor',decim, ...
    'SamplesPerFrame',samplesPerFrame, ...
    'OutputDataType','double', ...
    'EnableBurstMode',true, ...
    'NumFramesInBurst',1);

%% ================= 变量初始化 =================

CNT_END = 200;

count = 1;

Ber = [];

SNR_EVM_dB = [];

%% ================= 主循环 =================

while 1

    [rxSig,len] = step(radio);

    if len == 0
        continue
    end

    rxmimo2x2 = rxSig;

    figure_display = figure(1);

    subplot(2,3,1)

    set(figure_display,'units','normalized','position',[0.2473 0.1257 0.6258 0.5736]);

    plot(abs(rxmimo2x2));

    xlabel('采样点','FontWeight','bold','FontName','fangsong','FontSize',15);
    ylabel('幅度','FontWeight','bold','FontName','fangsong','FontSize',15);
    title('OFDM回波信号','FontWeight','bold','FontName','fangsong','FontSize',15);

%% ================= PSS同步 =================

    pss_vector = 1:length(PSS_code2);

    auto_correlation_vector = zeros(length(rxmimo2x2),1);

    for delay = 0:length(rxmimo2x2)-length(pss_vector)

        window = delay + pss_vector;

        auto_correlation_vector(delay+1) = ...
            PSS_transpose.' * (rxmimo2x2(window));

    end

    subplot(2,3,2)

    plot(abs(auto_correlation_vector))

    xlabel('采样点','FontWeight','bold','FontName','fangsong','FontSize',15);
    ylabel('相关性','FontWeight','bold','FontName','fangsong','FontSize',15);
    title('前导码同步','FontWeight','bold','FontName','fangsong','FontSize',15);

    [~,p1] = max(abs(auto_correlation_vector));

    temp = auto_correlation_vector;

    temp(p1) = 0;

    [~,p2] = max(abs(temp));

    p = min(p1,p2);

%% ================= OFDM符号定位 =================

    pos_cp1 = p + length(PSS_code2) - 1;

    pos_cp = pos_cp1 : N_symbol : ...
        (pos_cp1 + N_symbol * (Nd * spread_ratio - 1));

    if pos_cp(Nd * spread_ratio) + N_cp + N_fft - 1 > length(rxmimo2x2)
        continue
    end

%% ================= 去CP =================

    dcp = zeros(N_fft, data_col * spread_ratio);

    for i = 1:data_col * spread_ratio

        dcp_temp = rxmimo2x2( ...
            pos_cp(i) + N_cp : ...
            pos_cp(i) + N_cp + N_fft - 1);

        dcp(:,i) = dcp_temp;

    end

%% ================= FFT =================

    fft_data = fft(dcp, N_fft);

    data3 = fft_data(pos_row : pos_row + data_row + pilot_num - 1, :);

%% ================= 信道估计 =================

    Rx_pilot = data3(P_f_station,:);

    h = Rx_pilot ./ pilot_seq;

    H = interp1(P_f_station', h, data_station', 'linear','extrap');

%% ================= 均衡 =================

    data_aftereq = data3(data_station,:) ./ H;

    data_aftereq0 = data_aftereq;

%% ================= 并串 =================

    data_aftereq = reshape(data_aftereq,[],1);

%% ================= EVM SNR =================

    SNR_EVM_dB(count) = constellation_evm(data_aftereq.',2);

    mean_a = zeros(data_row,1);

    for i = 1:data_row
        mean_a(i) = mean(abs(data_aftereq0(i,:)));
    end

    subplot(2,3,3)

    plot(mean_a)

    xlabel('子载波索引','FontWeight','bold','FontName','fangsong','FontSize',15);
    ylabel('幅度','FontWeight','bold','FontName','fangsong','FontSize',15);
    title('频域接收信号(数据部分）','FontWeight','bold','FontName','fangsong','FontSize',15);

%% ================= 解扩 =================

    demspread_data = despread(data_aftereq0, code);

    demspread_data_row = reshape(demspread_data, [], 1);

%% ================= QPSK解调 =================

    demodulation_data = pskdemod(demspread_data_row, M, pi/M);

    De_data1 = reshape(demodulation_data, [], 1);

    De_data2 = de2bi(De_data1);

    De_Bit = reshape(De_data2', 1, []);

%% ================= 维特比译码 =================

    rx_c_de = vitdec(De_Bit, trellis, tblen, 'trunc', 'hard');

%% ================= BER =================

    [err(count), Ber(count)] = biterr(rx_c_de, P_data);

    fprintf('Frame %d  BER = %f\n', count, Ber(count));

%% ================= SNR显示 =================

    subplot(2,3,4)

    plot(SNR_EVM_dB)

    str3 = ['SNR_{E}：', num2str(SNR_EVM_dB(count))];

    text((1+count)/2,6,str3)

    xlabel('采集次数','FontWeight','bold','FontName','fangsong','FontSize',15);
    ylabel('信噪比（dB）','FontWeight','bold','FontName','fangsong','FontSize',15);
    title('直接计算接收信噪比','FontWeight','bold','FontName','fangsong','FontSize',15);

%% ================= BER显示 =================

    subplot(2,3,6)

    plot(Ber)

    hold on

    str = ['误码率：', num2str(Ber(count))];

    text((1+count)/2,0.3,str)

    ylim([0 0.5])

    hold off

    xlabel('采集次数','FontWeight','bold','FontName','fangsong','FontSize',15);
    ylabel('误码率','FontWeight','bold','FontName','fangsong','FontSize',15);
    title('接收信号误码率','FontWeight','bold','FontName','fangsong','FontSize',15);

%% ================= 循环控制 =================

    count = count + 1;

    if count > CNT_END

        fprintf('\n平均BER = %f\n', mean(Ber));

        return

    end

end
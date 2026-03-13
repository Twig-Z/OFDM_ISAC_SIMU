function SystemSetup_simu
%% ================================================================
%  SystemSetup_simu.m
%  功能：OFDM系统统一参数配置文件
%  被 OFDM_corr_rx.m 和 OFDM_corr_tx.m 共同调用
%
%  参数确定顺序（不可颠倒）：
%  masterClock → decim/interp → fs → N_fft → delta_f → N_sc → B(结果)
%  ================================================================

%% ================= 依赖文件检查 =================

[fList1, ~] = matlab.codetools.requiredFilesAndProducts('OFDM_proc_tx.m');
[fList2, ~] = matlab.codetools.requiredFilesAndProducts('OFDM_proc_rx.m');
fList = union(fList1, fList2);
disp('USRP收发调用的非系统函数/文件有：');
disp(fList')

%% ================= USRP硬件参数 =================
% 注意：masterClock和fs由硬件决定，不可随意更改

fc          = 915e6;        % 载波频率，可调范围：70MHz～6GHz
gain_tx     = 70;           % 发射增益，可调
gain_rx     = 30;           % 接收增益，可调

samplesPerFrame = 64000;    % 每帧采样点数

masterClock = 61.44e6;      % 硬件基准时钟，不可更改，硬件决定

decim       = 4;            % 抽取因子（RX），建议取2的幂次
interp      = decim;            % 插值因子（TX），建议取2的幂次，需与decim一致

fs = masterClock / decim;   % 实际采样率，由硬件决定，不可单独更改
                            % fs = 61.44MHz / 8 = 7.68MHz

%% ================= OFDM基础参数 =================
% 注意：delta_f由fs和N_fft共同决定，N_fft必须是2的幂次
mu      = 0;                        % 子载波间隔缩放因子（5G NR定义）
delta_f = 15000 * 2^mu;             % 子载波间隔15kHz，LTE/5G NR标准，建议不改

N_fft = fs / delta_f;               % FFT点数 = 7.68MHz / 15000 = 512
assert(N_fft == floor(N_fft),       'N_fft不是整数，检查fs与delta_f是否匹配');
assert(bitand(N_fft, N_fft-1) == 0, 'N_fft必须是2的幂次');

N_cp     = floor(0.07 * N_fft);     % 循环前缀长度 = 35
N_symbol = N_fft + N_cp;            % 单个OFDM符号长度 = 547

%% ================= 带宽与子载波参数 =================
% 注意：B是计算结果，不是输入量

eta          = 0.4;                         % 频谱效率因子，留出滚降保护
B_leackage   = 0.1 * fs;                    % 保护带宽 = 0.768MHz（fs的10%）

stage        = 1;                           % m序列扩频级数
spread_ratio = 2^stage - 1;                 % 扩频比 = 1（stage=1时无扩频）

B0   = fs - B_leackage;                     % 单路可用带宽 = 6.912MHz
B    = B0 * spread_ratio + B_leackage;      % 扩频后系统总带宽（结果量，勿作输入）

N_sc0 = floor(B0 / delta_f);               % 不含eta的子载波估计 = 460
assert(N_fft >= N_sc0, ...
    sprintf('N_fft(%d)装不下子载波N_sc0(%d)，请减小B_leackage或增大decim', ...
    N_fft, N_sc0));

N_sc    = floor(eta * B0 / delta_f);        % 数据+导频子载波总数 = 368
pos_row = ceil(B_leackage / delta_f / 2);   % FFT中数据起始bin偏移 = 26

%% ================= 子载波结构参数 =================

P_f_inter = 3;                              % 导频间隔（每3个子载波插1个导频）

% 导频子载波位置
P_f_station = 1 : P_f_inter : N_sc;
pilot_num   = length(P_f_station);          % 导频数量

% 数据子载波位置（非导频位置）
data_station = [];
for img = 1:N_sc
    if mod(img, P_f_inter) ~= 1
        data_station = [data_station, img];
    end
end
data_row = length(data_station);            % 数据子载波数量

P_f = sqrt(1/2) * (1 + 1i);                % 导频符号（QPSK点）

%% ================= 信道编码与调制参数 =================

M                = 4;           % 调制阶数：QPSK
R_channelcoding  = 1/2;         % 卷积码码率
L                = 7;           % 卷积码约束长度
tblen            = 6 * L;       % 维特比译码回溯长度 = 42
Nd               = 14;          % 每次传输的OFDM符号数
data_col = Nd;
% 每帧数据比特数
N_RE_data = data_row * data_col;
N_b_data  = N_RE_data * log2(M) * R_channelcoding;

%% ================= 扩频参数 =================

ptap1 = [1];                    % m序列生成多项式tap
regi1 = [1];                    % m序列初始寄存器状态

%% ================= PSS同步序列参数 =================

PSS_stage = 9;                          % PSS用m序列级数
PSS_ptap1 = [4 9];                      % PSS生成多项式tap
PSS_regi1 = [1 1 0 1 1 0 1 1 1];       % PSS初始寄存器状态

%% ================= 其他参数 =================

N_zeros = 10000;                % 发送前导零样本数（用于USRP稳定）

%% ================= 参数汇总打印 =================

fprintf('\n========== 系统参数汇总 ==========\n');
fprintf('硬件层:  masterClock=%.2fMHz, decim/interp=%d\n', masterClock/1e6, decim);
fprintf('采样层:  fs=%.3fMHz\n', fs/1e6);
fprintf('OFDM层:  delta_f=%.0fHz, N_fft=%d, N_cp=%d\n', delta_f, N_fft, N_cp);
fprintf('带宽层:  B0=%.3fMHz, B=%.3fMHz\n', B0/1e6, B/1e6);
fprintf('子载波:  N_sc=%d, pilot=%d, data=%d, pos_row=%d\n', N_sc, pilot_num, data_row, pos_row);
fprintf('编码层:  M=%d, R=%.1f, L=%d, Nd=%d\n', M, R_channelcoding, L, Nd);
fprintf('扩频层:  stage=%d, spread_ratio=%d\n', stage, spread_ratio);
fprintf('===================================\n\n');

% 获取当前程序所在路径
currentFolder = fileparts(mfilename('fullpath'));
save(fullfile(currentFolder, 'parameter_ttt.mat'));  

end
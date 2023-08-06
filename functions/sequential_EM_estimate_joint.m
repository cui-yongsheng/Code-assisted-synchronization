function [Est_phase,Est_freq,Est_time,LLR_rec] = sequential_EM_estimate_joint(QPSK_frame_sample,symbol_noise_var,H,sps)
%EM_ESTIMATE 基于EM算法LDPC估计载波频率(采用序贯策略)
%   输入：        1.           QPSK_frame_sample              QPSK符号帧(sps采样）
%                 2.           symbol_noise_var        噪声方差
%                 3.           H                       LDPC校验矩阵
%                 4.           sps                     采样率
%
%   输出：        1.           Est_phase   相位估计值(以向量形式存储估计值迭代变化)
%                 2.           Est_freq    频率估计值(以向量形式存储估计值迭代变化)
%                 3.           Est_time    定时估计值(以向量形式存储估计值迭代变化)
%                 4.           LLR_rec     译码后软信息

%% LDPC译码器辅助参数配置
H_info=gen_H_info(H);
H_weight = H_info.H_weight;
RowIndex = H_info.RowIndex;
ColWeight = H_info.ColWeight;
N = length(ColWeight);                                % 校验矩阵列数 

%% 内置参数
EM_iterNum = 20;
learning_rate_index = (0:EM_iterNum-1)./EM_iterNum.*pi;
learning_rate = (cos(learning_rate_index)+1).*1e-4;         %更新步长
gama = 1e-4;                                                %近似位置
f1 = -4e-5;
f2 = 4e-5;

%% QPSK软解调
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK解调器(软解调，并根据输入的噪声功率，
    'DecisionMethod','Approximate log-likelihood ratio', ...    % 在假定复信号功率为1的条件下进行LLR归一化
    'VarianceSource','Input port');

%% 变量节点LLR初始化
QPSK_frame = QPSK_frame_sample(1:sps:end);
y = qpskDemod(QPSK_frame,symbol_noise_var);
v = zeros(H_weight+1,1);                                % 变量节点LLR向量
col_index = 1;
for col=1:N
    v(col_index:col_index+ColWeight(col)-1) = y(col);   % 基于信道先验LLR初始化变量节点LLR向量
    col_index = col_index+ColWeight(col);               % 列索引更新
end

%% EM循环迭代
Est_freq = zeros(1,EM_iterNum+1);
Est_phase = zeros(1,EM_iterNum+1);
Est_time = zeros(1,EM_iterNum+1);
for i=1:EM_iterNum
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',-1*Est_freq(i),'PhaseOffset',-1*Est_phase(i),'SampleRate',sps);
    QPSK_frame_trail_sample=ErrIterator(QPSK_frame_sample);
    QPSK_frame_trail=TEmodel(QPSK_frame_trail_sample,Est_time(i));
    QPSK_frame_trail=QPSK_frame_trail(1:sps:end);
    LLR=qpskDemod(QPSK_frame_trail,symbol_noise_var);
    [LLR_rec,~,v] = fun_BP_decode_mex(H_weight,RowIndex,ColWeight,LLR,1,v);
    BPSK_rec=tanh(LLR_rec/2);
    QPSK_rec=1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end);                        %软符号重构
    [F_max,Angle_max]=give_CZT_graph(1,QPSK_frame.*conj(QPSK_rec),f1,f2);   
    Est_freq(i+1) = F_max;
    Est_phase(i+1) = Angle_max/pi*180;
    QPSK_frame_trail_right=TEmodel(QPSK_frame_trail_sample,Est_time(i)+gama);
    QPSK_frame_trail_right=QPSK_frame_trail_right(1:sps:end);
    QPSK_frame_trail_left=TEmodel(QPSK_frame_trail_sample,Est_time(i)-gama);
    QPSK_frame_trail_left=QPSK_frame_trail_left(1:sps:end);
    QPSK_frame_trail_gradient=(QPSK_frame_trail_right-QPSK_frame_trail_left)/gama;
    s_curve = real(sum(conj(QPSK_rec).*QPSK_frame_trail_gradient));
    Est_time(i+1) = Est_time(i)+learning_rate(i)*s_curve*sps;      %唐发建
end
Est_time=-Est_time;
end


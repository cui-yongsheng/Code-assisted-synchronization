function [phasePre,freqPre,timePre] = phaseFreqTimeSearch(QPSK_frame_sample,symbol_noise_var,H,sps)
%phaseFreqTimeSearch 通过多点预搜索扩大收敛范围
%   输入：        1.           QPSK_frame_sample       QPSK符号帧(sps采样)
%                 2.           symbol_noise_var        噪声方差
%                 3.           H                       LDPC校验矩阵
%                 4.           sps                     采样率
%
%
%   输出：        1.           phasePre            相位搜索值
%                 2.           freqPre             频率搜索值
%                 3.           timePre             时延搜索值

%% Todo: 校验子范数计算方法优化

%% 内置参数     
BP_times=1;                     % BP译码迭代次数
Threshold=0.1;                  % 阈值(判断校验子范数是否有效)
list_length=5;                  % 候选列表元素个数

%% 相频预搜索参数设置
G_phi=52;                       % 相位搜索间隔                                     
N_phi=3;                        % 相位搜索点数
G_f=1e-4;                       % 频率搜索间隔
N_f=4;                          % 频率搜索点数
G_t=sps/4;                      % 定时搜索间隔
N_t=2;                          % 定时搜索点数

%% 通信工具箱对象实例化
ldpcDec_for_Fun = comm.LDPCDecoder('ParityCheckMatrix',H, ...
    'DecisionMethod','Soft decision', ...
    'MaximumIterationCount',BP_times, ...
    'OutputValue','Whole codeword','FinalParityChecksOutputPort',true);     % LDPC译码器
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK解调器(软解调，并根据输入的噪声功率，
    'DecisionMethod','Approximate log-likelihood ratio', ...    % 在假定复信号功率为1的条件下进行LLR归一化
    'VarianceSource','Input port');

%% 获取所有选择点的目标函数值
access_phase_vec=(-N_phi:N_phi)*G_phi;          % 相位查找列表
access_freq_vec=(-N_f:N_f)*G_f;                 % 频率查找列表
access_time_vec=(-N_t:N_t)*G_t;                 % 定时查找列表
Object_Fun_CMF=zeros(length(access_phase_vec),length(access_freq_vec),length(access_time_vec));        %校验子范数存储空间
for timeErrList_index=1:length(access_time_vec)
    for freqErrList_index=1:length(access_freq_vec)
        for phaseErrList_index=1:length(access_phase_vec)
            CMF_ErrIterator=comm.PhaseFrequencyOffset('PhaseOffset',access_phase_vec(phaseErrList_index),'FrequencyOffset',access_freq_vec(freqErrList_index),'SampleRate',sps);
            rxFrame_checked_sample=CMF_ErrIterator(QPSK_frame_sample);
            rxFrame_checked_sample=TEmodel(rxFrame_checked_sample,access_time_vec(timeErrList_index));
            rxFrame_checked=rxFrame_checked_sample(1:sps:end);
            rxLLR_bit_stream=qpskDemod(rxFrame_checked,symbol_noise_var);
            Hard_vector=(rxLLR_bit_stream<=0);
            % 校验子范数
            Object_Fun_CMF(phaseErrList_index,freqErrList_index,timeErrList_index)=sum(1-double(rem(sum(H(:,Hard_vector),2),2)));   
        end
    end
end
%% 对所有点比较以给出初步搜索结果
[~,I_z]=max(max(max(Object_Fun_CMF)));
[~,I_y]=max(max(Object_Fun_CMF(:,:,I_z)));
[~,I_x]=max(Object_Fun_CMF(:,I_y,I_z));
phasePre=access_phase_vec(I_x);
freqPre=access_freq_vec(I_y);
timePre=access_time_vec(I_z);

%% 判断初步结果是否符合要求，符号要求直接返回，否则计算目标函数
search_list = sort(max(max(Object_Fun_CMF)));
wait_list = zeros(1,list_length);                             % 候选列表
if (search_list(1)-search_list(2))/search_list(1) < Threshold           % 条件成立为不符合要求
    % 返回指定元素下标
    [~,sorted_indices] = maxk(Object_Fun_CMF(:),5);                     % 校验子范数最大值下标
    [row,col,depth] = ind2sub(size(Object_Fun_CMF),sorted_indices);
    % 对候选列表中元素重新计算目标函数值
    for index=1:list_length
        CMF_ErrIterator=comm.PhaseFrequencyOffset('PhaseOffset',access_phase_vec(row(index)),'FrequencyOffset',access_freq_vec(col(index)),'SampleRate',sps);
        rxFrame_checked_sample=CMF_ErrIterator(QPSK_frame_sample);
        rxFrame_checked_sample=TEmodel(rxFrame_checked_sample,access_time_vec(depth(index)));
        rxFrame_checked=rxFrame_checked_sample(1:sps:end);
        rxLLR_bit_stream=qpskDemod(rxFrame_checked,symbol_noise_var);
        % 计算目标函数
        LLR_rec=ldpcDec_for_Fun(rxLLR_bit_stream);
        wait_list(index)=0.5*mean(log(1+cosh(LLR_rec)));
    end
    [~,index] = max(wait_list);
    phasePre=access_phase_vec(row(index));
    freqPre=access_freq_vec(col(index));
    timePre=access_time_vec(depth(index));
end
end
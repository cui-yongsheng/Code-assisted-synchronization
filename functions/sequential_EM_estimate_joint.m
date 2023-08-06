function [Est_phase,Est_freq,Est_time,LLR_rec] = sequential_EM_estimate_joint(QPSK_frame_sample,symbol_noise_var,H,sps)
%EM_ESTIMATE ����EM�㷨LDPC�����ز�Ƶ��(����������)
%   ���룺        1.           QPSK_frame_sample              QPSK����֡(sps������
%                 2.           symbol_noise_var        ��������
%                 3.           H                       LDPCУ�����
%                 4.           sps                     ������
%
%   �����        1.           Est_phase   ��λ����ֵ(��������ʽ�洢����ֵ�����仯)
%                 2.           Est_freq    Ƶ�ʹ���ֵ(��������ʽ�洢����ֵ�����仯)
%                 3.           Est_time    ��ʱ����ֵ(��������ʽ�洢����ֵ�����仯)
%                 4.           LLR_rec     ���������Ϣ

%% LDPC������������������
H_info=gen_H_info(H);
H_weight = H_info.H_weight;
RowIndex = H_info.RowIndex;
ColWeight = H_info.ColWeight;
N = length(ColWeight);                                % У��������� 

%% ���ò���
EM_iterNum = 20;
learning_rate_index = (0:EM_iterNum-1)./EM_iterNum.*pi;
learning_rate = (cos(learning_rate_index)+1).*1e-4;         %���²���
gama = 1e-4;                                                %����λ��
f1 = -4e-5;
f2 = 4e-5;

%% QPSK����
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK�����(������������������������ʣ�
    'DecisionMethod','Approximate log-likelihood ratio', ...    % �ڼٶ����źŹ���Ϊ1�������½���LLR��һ��
    'VarianceSource','Input port');

%% �����ڵ�LLR��ʼ��
QPSK_frame = QPSK_frame_sample(1:sps:end);
y = qpskDemod(QPSK_frame,symbol_noise_var);
v = zeros(H_weight+1,1);                                % �����ڵ�LLR����
col_index = 1;
for col=1:N
    v(col_index:col_index+ColWeight(col)-1) = y(col);   % �����ŵ�����LLR��ʼ�������ڵ�LLR����
    col_index = col_index+ColWeight(col);               % ����������
end

%% EMѭ������
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
    QPSK_rec=1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end);                        %������ع�
    [F_max,Angle_max]=give_CZT_graph(1,QPSK_frame.*conj(QPSK_rec),f1,f2);   
    Est_freq(i+1) = F_max;
    Est_phase(i+1) = Angle_max/pi*180;
    QPSK_frame_trail_right=TEmodel(QPSK_frame_trail_sample,Est_time(i)+gama);
    QPSK_frame_trail_right=QPSK_frame_trail_right(1:sps:end);
    QPSK_frame_trail_left=TEmodel(QPSK_frame_trail_sample,Est_time(i)-gama);
    QPSK_frame_trail_left=QPSK_frame_trail_left(1:sps:end);
    QPSK_frame_trail_gradient=(QPSK_frame_trail_right-QPSK_frame_trail_left)/gama;
    s_curve = real(sum(conj(QPSK_rec).*QPSK_frame_trail_gradient));
    Est_time(i+1) = Est_time(i)+learning_rate(i)*s_curve*sps;      %�Ʒ���
end
Est_time=-Est_time;
end


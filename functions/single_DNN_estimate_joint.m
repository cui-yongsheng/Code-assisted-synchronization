function [Est_phase,Est_freq,Est_time] = single_DNN_estimate_joint(QPSK_frame_sample,symbol_noise_var,H,sps)
%DNN_ESTIMATE ����DNN�����ز�Ƶ��
%   ���룺        1.           QPSK_frame              QPSK����֡(sps����)
%                 2.           symbol_noise_var        ��������
%                 3.           H                       LDPCУ�����
%                 4.           sps                     ������
%
%   �����        1.           Est_phase   ��λ����ֵ(��������ʽ�洢����ֵ�����仯)
%

%% ���ò���
DNN_iterNum = 1;            % DNN��������
BP_Num = 2;

%% LDPC������
ldpcDec_for_Fun = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Soft decision','MaximumIterationCount',BP_Num,'OutputValue','Whole codeword');

%% QPSK����
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK�����(������������������������ʣ�
    'DecisionMethod','Approximate log-likelihood ratio', ...    % �ڼٶ����źŹ���Ϊ1�������½���LLR��һ��
    'VarianceSource','Input port');

%% DNNѭ������
Est_phase = zeros(1,DNN_iterNum+1);
Est_freq = zeros(1,DNN_iterNum+1);
Est_time = zeros(1,DNN_iterNum+1);
predict_phase = zeros(1,DNN_iterNum+1);
predict_freq = zeros(1,DNN_iterNum+1);
predict_time = zeros(1,DNN_iterNum+1);
for i=1:DNN_iterNum
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',-1*Est_freq(i),....
    'PhaseOffset',-1*Est_phase(i),'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(QPSK_frame_sample);
    sample_time = 1:length(rxSig_checked_sample);
    sample_time_adjust=1+Est_time(i):length(rxSig_checked_sample)+Est_time(i);
    rxSig_checked_sample=(interp1(sample_time,rxSig_checked_sample',sample_time_adjust,"linear"))';    %����������ֵ
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);
    
    %% ����Ŀ�꺯��0
    LLR_vector_self=ldpcDec_for_Fun(LLR);
    Fun_d0=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d1_p=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d1_f=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    rxSig_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxSig_gradient_sample=(rxSig_adjust_sample-rxSig_checked_sample)./0.01;
    rxSig_gradient_trail=rxSig_gradient_sample(1:sps:end);
    Fun_d1_t=real(sum(rxSig_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�
    
    %% �ն˴���1 ������λ�ã�
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',-1*Est_freq(i)+1e-5,....
    'PhaseOffset',-1*Est_phase(i)+5,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(QPSK_frame_sample);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,Est_time(i)+1);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);
    
    %% ����Ŀ�꺯��1
    LLR_vector_self=ldpcDec_for_Fun(LLR);
    Fun_d0_n=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d1_p_n=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d1_f_n=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    rxSig_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxSig_gradient_sample=(rxSig_adjust_sample-rxSig_checked_sample)./0.01;
    rxSig_gradient_trail=rxSig_gradient_sample(1:sps:end);
    Fun_d1_t_n=real(sum(rxSig_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�

    %% ������������
    X_data=[Fun_d0 Fun_d1_p Fun_d1_f Fun_d1_t Fun_d0_n Fun_d1_p_n Fun_d1_f_n Fun_d1_t_n];
    predict_phase(i) = NNF_p(X_data);
    predict_freq(i) = NNF_f(X_data)/1e6;
    predict_time(i) = NNF_t(X_data);
    Est_phase(i+1) = Est_phase(i)+predict_phase(i);
    Est_freq(i+1) = Est_freq(i)+predict_freq(i);
    Est_time(i+1) = Est_time(i)+predict_time(i);
    
end

end


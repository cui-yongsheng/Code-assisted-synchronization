function txSig = TEmodel(txSig,timingErr)
%TEmodel 添加定时偏移。要求输入为列向量,偏移样本点数应小于sps/2。
%   输入：
%               1.txSig             发送信号
%               2.timingErr         定时偏差样本点数
%
%   输出：
%               1.txSig             输出信号
%   示例：
%   输入：
%               txSig=0 1 2 3 4 5
%               timingErr=1
%   输出：
%               txSig=1 2 3 4 5 6
%  

%% Todo: 添加频率时间偏移

%% 进行插值偏移
sample_time = (1:length(txSig))';
sample_time_adjust = sample_time+timingErr;
txSig = interp1(sample_time,txSig,sample_time_adjust,'spline');
end

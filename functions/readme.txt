1.NFF_f:        用于预测频率的DNN网络
2.NFF_p:        用于预测相位的DNN网络
3.NFF_t:        用于预测时延的DNN网络
4.phaseFreqTimeSearch:      预搜索扩大收敛范围
5.sequential_EM_estimate_joint：     预测频率、相位、时延，LDPC译码采用序贯策略 
6.single_DNN_estimate_joint：        基于DNN预测  
7.TEmodel:      使用插值实现的时延操作


TODO：
4.phaseFreqTimeSearch:      校验子范数计算方法
7.TEmodel:                  Todo: 添加频率时间偏移
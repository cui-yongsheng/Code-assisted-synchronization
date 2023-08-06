1.data.mat          10000帧，包含22个特征，3个标签。具体见sample_gen.m。译码迭代次数为6
2.sample_gen.m      生成数据
3.train.m           进行训练，并生成DNN函数。（数据输入部分需要调整，目前可方便调试）



不同点查找效果：

样本:5000
点数:4
BP:2

特征:ponit0 ponit1 ponit2 ponit3 look_fun_table look_fun_table1
效果:phase 11.7555    freq 30.4942    time 0.0033

特征:ponit0 ponit3 look_fun_table look_fun_table1
效果:phase 11.6372    freq 38.4208    time 0.0037

特征:ponit0 ponit3 
效果:phase 11.5643    freq 43.9362    time 0.0038

特征:ponit0 
效果:phase 12.4854    freq 48.5445    time 0.0061

样本:5000
点数:4
BP:4

特征:ponit0 
效果:phase 8.8530    freq 29.1212   time 0.0049

特征:ponit0 ponit3 
效果:phase 7.7748    freq 26.8700   time 0.0032

特征:ponit0 ponit3 look_fun_table look_fun_table1
效果:phase 7.5126    freq 24.3209    time 0.0032

样本:5000
点数:4
BP:6

特征:ponit0 ponit3 look_fun_table look_fun_table1
效果:phase 6.3762    freq 20.9167    time 0.0032

特征:ponit0 ponit3
效果:phase 6.7843    freq 19.0188    time 0.0031

特征:ponit0 
效果:phase 7.0782    freq 22.1655    time 0.0042

##################################################
选择特征
样本:5000
点数:4
BP:6
特征:ponit0 ponit3
效果:phase 6.7843    freq 19.0188    time 0.0031
##################################################

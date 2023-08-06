function [] = send_message(head,text)
%SEND_MESSAGE 此函数使用server酱向微信推送信息
%   输入：
%           1.head      标题
%           2.text      文件内容
%   输出：
%           

KEY = 'SCT206826TP6Y3Ya9rGWmEBmDpazSEXZjU';

if nargin==0
   head = 'matlab_program_over';
   text = '啦啦啦~';
elseif nargin==1
   text = '啦啦啦~';
end
URL = ['http://sc.ftqq.com/',KEY,'.send?text=',head,'&desp=',text];
try
    webread(URL);
catch 
    disp('次数超出限制');
end
end


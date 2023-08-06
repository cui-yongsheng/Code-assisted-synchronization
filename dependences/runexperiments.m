files = dir('*');
matchedFiles = {};
tic
for i = 1:length(files)
    filename = files(i).name;
    if startsWith(filename,'test') || endsWith(filename,'test.m')
        tic
        matchedFiles{end+1} = filename;
        disp(['当前文件:',filename])
        try
            run(filename)
        catch exception
            disp(['出现错误:',exception.name])
        end
        toc
    end
end
toc


log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun0'};

for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(8)
        fid = fopen([log_path 'MBA_' func_name '_' num2str(N) '.log'],'a+');
        for EPS = [5 7]
            runCbs(N, func_name, EPS, fid);
        end
        fclose(fid);
    end
end

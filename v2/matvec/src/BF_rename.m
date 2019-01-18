function y = BF_rename(Factor)
% This code renames the matrix factors of a butterfly matrix
%
% Copyright 2018 by Haizhao Yang

switch numel(fieldnames(Factor))
    case 5
        len = length(Factor.HTol);
        n = 2*(len+1)+1;
        y = cell(1,n);
        y{n} = Factor.V;
        
        for i=len:-1:1
            y{len+2+i} = Factor.HTol{i};
        end
        
        y{len+2} = Factor.M;
        
        for i=1:len
            y{len+2-i} = Factor.GTol{i};
        end
        
        y{1} = Factor.U;
    case 4
        len = length(Factor.U);
        n = 2*len+1;
        y = cell(1,n);
        
        for i=1:len
            y{2*len+2-i} = Factor.V{i};
        end
        
        y{len+1} = Factor.S;
        
        for i=len:-1:1
            y{i} = Factor.U{i};
        end
    case 3
        len = length(Factor.U);
        n = 2*len+1;
        y = cell(1,n);
        
        for i=1:len
            y{2*len+2-i} = Factor.V{i};
        end
        
        y{len+1} = Factor.S;
        
        for i=len:-1:1
            y{i} = Factor.U{i};
        end
end

end

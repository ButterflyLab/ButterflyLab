raw_path = '../test/log/';

func_name = 'fun0';

NGlist = [4, 6, 8];

outfid = cell(100,1);
for NG = NGlist
    outfid{NG} = fopen(['Table_' func_name '_1D_' num2str(NG) '.log'],...
        'w');
end

for N = 2.^(8:2:20)
    fid = fopen([raw_path 'Factor_' func_name '_1D_' num2str(N) '.log'],...
        'r');
    eline = fgets(fid); %----------------------------------------------
    while ischar(eline)
        tline   = fgets(fid);
        tN      = textscan(tline,'N                 : %d');
        tline   = fgets(fid);
        NG      = textscan(tline,'Chebyshev pts     : %d');
        tline   = fgets(fid);
        tol     = textscan(tline,'Tolerance         : %f');
        tline   = fgets(fid);
        relerr2 = textscan(tline,'Relative Error_2  : %f');
        tline   = fgets(fid);
        Rcomp   = textscan(tline,'Compression Ratio : %f');
        tline   = fgets(fid);
        Td      = textscan(tline,'Direct Time       : %f s');
        tline   = fgets(fid);
        Tr      = textscan(tline,'Running Time      : %f mins');
        tline   = fgets(fid);
        Tfact   = textscan(tline,'Factorization Time: %f mins');
        tline   = fgets(fid);
        Ta      = textscan(tline,'Applying Time     : %f s');
        fgets(fid); %----------------------------------------------
        tline = fgets(fid);
        eline = fgets(fid);
        fprintf(outfid{NG{1}},...
            '%7d,%2d & %.2e & %.2e & %.2e & %.2e & %.2e \\\\\n',...
            tN{1},NG{1},relerr2{1},Tfact{1},Td{1},Ta{1},Td{1}/Ta{1});
    end
    fclose(fid);
end

for NG = NGlist
    fclose(outfid{NG});
end

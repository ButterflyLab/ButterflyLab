function kidx = bf_prep(kk,kbox,npk1,npk2)

kidx = cell(npk1,npk2);
k1offset = kbox(1,1);
k1len = kbox(1,2)-kbox(1,1);
k2offset = kbox(2,1);
k2len = kbox(2,2)-kbox(2,1);
for i=1:npk1
    for j=1:npk2
        k1s = k1offset + (i-1)*k1len/npk1;
        k1e = k1offset + i*k1len/npk1;
        k2s = k2offset + (j-1)*k2len/npk2;
        k2e = k2offset + j*k2len/npk2;
        kidx{i,j} = find( kk(:,1)>=k1s & kk(:,1)<k1e ...
            & kk(:,2)>=k2s & kk(:,2)<k2e );
    end
end

end

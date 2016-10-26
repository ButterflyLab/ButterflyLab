function kidx = bf_prep(kk,kbox,npk)

kidx = cell(npk,1);
koffset = kbox(1);
klen = kbox(2)-kbox(1);
for i=1:npk
    ks = koffset + (i-1)*klen/npk;
    ke = koffset + i*klen/npk;
    kidx{i} = find( kk>=ks & kk<ke );
end

end
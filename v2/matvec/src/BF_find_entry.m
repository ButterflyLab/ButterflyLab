function k = BF_find_entry(a,b)
% find the positions of the elements of b in a
% TODO: use bisection method to speed up
k = zeros(size(b));
for cnt = 1:numel(b)
    if abs(cnt-1)<1e-10
        st = 1;
    else
        st = k(cnt-1)+1;
    end
    for cntc = st:numel(a)
        if abs(b(cnt)-a(cntc))<1e-10
            k(cnt) = cntc;
            break;
        end
    end
end


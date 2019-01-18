function [U,V] = BF_Phase_Correction_Mat_2D(U,V,posu,posv,disPosu,disPosv,disPosSubu,disPosSubv,tau)
% U = mod(Phi(:,posu),tau) and V = mod(Phi(posv,:)',tau)
% disPosu(:,1) contains the discontinuous positions along columns,
% including fake and true discontinuity
% disPosv(:,1) contains the discontinuous posisions along rows,
% including fake and true discontinuity
% disPosu(:,2) contains the type of discontinuous positions along columns,
% 0 (inner) and 1 (outer) for fake, 2 for true
% disPosv(:,2) contains the thpe of discontinuous posisions along rows,
% 0 (inner) and 1 (outer) for fake, 2 for true
% posu(disPosSubu) are the the positons for true discontinuity along
% columns
% posv(disPosSubv) are the the positons for true discontinuity along rows
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.


% U for columns and V for rows, all talk skinny matrices
[mu,nu] = size(U);
[mv,nv] = size(V);

% correct first row in each continuous sector of rows
for cntd = 1:numel(disPosSubu)
    posc = disPosSubu(cntd);
    V(:,posc) = BF_Phase_Correction_Vec_1_2D(V(:,posc),disPosv,tau);
end

% correct first three columns in each continuous sector of columns
for cntd = 1:numel(disPosv)
    for cnt1 = 1:3
        pos = cnt1+disPosSubv(cntd)-1;
        val = V(disPosv(cntd)+cnt1-1,disPosSubu);
        U(:,pos) = BF_Phase_Correction_Vec_2(U(:,pos),tau,disPosu,val);
    end
    % for loop for each continuous sector of rows
    for cntdd = 1:numel(disPosu)
        app = U(disPosu(cntdd),disPosSubv(cntd)+1) + U(disPosu(cntdd)+1,disPosSubv(cntd)) - U(disPosu(cntdd),disPosSubv(cntd));
        slope = round((app-U(disPosu(cntdd)+1,disPosSubv(cntd)+1))/tau);
        if cntdd == numel(disPosu)
            len = mu-disPosu(cntdd)+1;
            ed = mu;
        else
            len = disPosu(cntdd+1)-disPosu(cntdd);%+1;
            ed = disPosu(cntdd+1)-1;
        end
        st = disPosu(cntdd);
        vecAdd = tau*(0:(len-1))*slope;
        U(st:ed,disPosSubv(cntd)+1) = U(st:ed,disPosSubv(cntd)+1) + vecAdd';
        app = U(disPosu(cntdd),disPosSubv(cntd)+2) + U(disPosu(cntdd)+1,disPosSubv(cntd)+1) - U(disPosu(cntdd),disPosSubv(cntd)+1);
        slope = round((app-U(disPosu(cntdd)+1,disPosSubv(cntd)+2))/tau);
        vecAdd = tau*(0:(len-1))*slope;
        U(st:ed,disPosSubv(cntd)+2) = U(st:ed,disPosSubv(cntd)+2) + vecAdd';
    end
end

%correct all other rows
for cntd = 1:numel(disPosu)
    if cntd == numel(disPosu)
        ed = nv;
    else
        ed = disPosSubu(cntd+1)-1;
    end
    for cnt1 = disPosSubu(cntd)+1:ed
        pos1 = disPosv; pos2 = posv(cnt1);
        val1 = U(pos2,disPosSubv); val2 = U(pos2,disPosSubv+1); val3 = U(pos2,disPosSubv+2);
        scnDer = (val1+val3-2*val2);
        V(:,cnt1) = BF_Phase_Correction_Vec_3(V(:,cnt1),tau,pos1,val1,pos1+1,val2,pos1+2,val3,scnDer);
    end
end

% correct all other columns
for cntd = 1:numel(disPosv)
    if cntd == numel(disPosv)
        ed = nu;
    else
        ed = disPosSubv(cntd+1)-1;
    end
    for cnt1 = disPosSubv(cntd)+3:ed
        pos1 = disPosu; pos2 = posu(cnt1);
        val1 = V(pos2,disPosSubu); val2 = V(pos2,disPosSubu+1); val3 = V(pos2,disPosSubu+2);
        scnDer = (val1+val3-2*val2);
        U(:,cnt1) = BF_Phase_Correction_Vec_3(U(:,cnt1),tau,pos1,val1,pos1+1,val2,pos1+2,val3,scnDer);
    end
end
end



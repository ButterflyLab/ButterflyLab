function D = deserialize(fid, type)
  
  switch type{1}
   case 'char'
    D = fread(fid, 1, 'char');
   case 'int'
    D = fread(fid, 1, 'int');
   case 'double'
    D = fread(fid, 1, 'double');
   case 'cpx'
    T = fread(fid, 2, 'double');
    D = T(1) + i*T(2);
   case 'Index2'
    D = fread(fid, 2, 'int');
   case 'Point2'
    D = fread(fid, 2, 'double');
   case 'Index3'
    D = fread(fid, 3, 'int');
   case 'Point3'
    D = fread(fid, 3, 'double');
   case 'vector'
    m = fread(fid, 1, 'int');
    D = cell(m,1);
    for k=1:numel(D)
      D{k} = deserialize(fid, type{2});
    end
   case 'map'
    m = fread(fid, 1, 'int');
    D = cell(m,2);
    for k=1:size(D,1)
      D{k,1} = deserialize(fid, type{2});
      D{k,2} = deserialize(fid, type{3});
    end
   case 'pair'
    D = cell(1,2);
    D{1} = deserialize(fid, type{2});
    D{2} = deserialize(fid, type{3});
   case 'tuple'
    D = cell(1,numel(type)-1);
    for k=1:numel(D)
      D{k} = deserialize(fid, type{k+1});
    end
   case 'BolNumVec'
    m = fread(fid, 1, 'int');
    D = fread(fid, m, 'char');
    D = reshape(D, [m,1]);
   case 'BolNumMat'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    D = fread(fid, m*n, 'char');
    D = reshape(D, [m,n]);
   case 'BolNumTns'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    p = fread(fid, 1, 'int');
    D = fread(fid, m*n*p, 'char');
    D = reshape(D, [m,n,p]);
   case 'IntNumVec'
    m = fread(fid, 1, 'int');
    D = fread(fid, m, 'int');
    D = reshape(D, [m,1]);
   case 'IntNumMat'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    D = fread(fid, m*n, 'int');
    D = reshape(D, [m,n]);
   case 'IntNumTns'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    p = fread(fid, 1, 'int');
    D = fread(fid, m*n*p, 'int');
    D = reshape(D, [m,n,p]);
   case 'DblNumVec'
    m = fread(fid, 1, 'int');
    D = fread(fid, m, 'double');
    D = reshape(D, [m,1]);
   case 'DblNumMat'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    D = fread(fid, m*n, 'double');
    D = reshape(D, [m,n]);
   case 'DblNumTns'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    p = fread(fid, 1, 'int');
    D = fread(fid, m*n*p, 'double');
    D = reshape(D, [m,n,p]);
   case 'CpxNumVec'
    m = fread(fid, 1, 'int');
    T = fread(fid, 2*m, 'double');
    rD = reshape(T(1:2:end), [m,1]);
    iD = reshape(T(2:2:end), [m,1]);
    D = rD + i*iD;
   case 'CpxNumMat'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    T = fread(fid, 2*m*n, 'double');
    rD = reshape(T(1:2:end), [m,n]);
    iD = reshape(T(2:2:end), [m,n]);
    D = rD + i*iD;
   case 'CpxNumTns'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    p = fread(fid, 1, 'int');
    T = fread(fid, 2*m*n*p, 'double');
    rD = reshape(T(1:2:end), [m,n,p]);
    iD = reshape(T(2:2:end), [m,n,p]);
    D = rD + i*iD;
   case 'NumVec'
    m = fread(fid, 1, 'int');
    D = cell(m,1);
    for k=1:numel(D)
      D{k} = deserialize(fid, type{2});
    end
   case 'NumMat'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    D = cell(m,n);
    for k=1:numel(D)
      D{k} = deserialize(fid, type{2});
    end
   case 'NumTns'
    m = fread(fid, 1, 'int');
    n = fread(fid, 1, 'int');
    p = fread(fid, 1, 'int');
    D = cell(m,n,p);
    for k=1:numel(D)
      D{k} = deserialize(fid, type{2});
    end
   case 'vectorPoint3'
    m = fread(fid,1,'int');
    D = fread(fid,3*m,'double');
    D = reshape(D, [3,m]);
   case 'vectorIndex3'
    m = fread(fid,1,'int');
    D = fread(fid,3*m,'int');
    D = reshape(D, [3,m]);
   otherwise
    error('wrong');
  end
  
  
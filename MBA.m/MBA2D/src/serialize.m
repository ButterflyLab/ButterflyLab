function serialize(fid, D, type)
  
  switch type{1}
   case 'char'
    fwrite(fid, D, 'char');
   case 'int'
    fwrite(fid, D, 'int');
   case 'double'
    fwrite(fid, real(D), 'double');
   case 'cpx'
    T = [real(D) imag(D)];
    fwrite(fid, real(T), 'double');
   case 'Index2'
    fwrite(fid, D, 'int');
   case 'Point2'
    fwrite(fid, D, 'double');
   case 'Index3'
    fwrite(fid, D, 'int');
   case 'Point3'
    fwrite(fid, D, 'double');
   case 'vector'
    fwrite(fid, size(D,1), 'int');
    for k=1:numel(D)
      serialize(fid, D{k}, type{2});
    end
   case 'map'
    fwrite(fid, size(D,1), 'int');
    for k=1:size(D,1)
      serialize(fid, D{k,1}, type{2});
      serialize(fid, D{k,2}, type{3});
    end
   case 'pair'
    serialize(fid, D{1}, type{2});
    serialize(fid, D{2}, type{3});
   case 'tuple'
    for k=1:numel(D)
      serialize(fid, D{k}, type{k+1});
    end
   case 'BolNumVec'
    fwrite(fid, numel(D), 'int');
    fwrite(fid, D, 'char');
   case 'BolNumMat'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, D, 'char');
   case 'BolNumTns'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, size(D,3), 'int');
    fwrite(fid, D, 'char');
   case 'IntNumVec'
    fwrite(fid, numel(D), 'int');
    fwrite(fid, D, 'int');
   case 'IntNumMat'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, D, 'int');
   case 'IntNumTns'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, size(D,3), 'int');
    fwrite(fid, D, 'int');
   case 'DblNumVec'
    fwrite(fid, numel(D), 'int');
    fwrite(fid, D, 'double');
   case 'DblNumMat'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, D, 'double');
   case 'DblNumTns'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, size(D,3), 'int');
    fwrite(fid, D, 'double');
   case 'CpxNumVec'
    fwrite(fid, numel(D), 'int');
    rD = real(D);    iD = imag(D);
    T = [rD(:)'; iD(:)'];
    fwrite(fid, T, 'double');
   case 'CpxNumMat'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    rD = real(D);    iD = imag(D);
    T = [rD(:)'; iD(:)'];
    fwrite(fid, T, 'double');
   case 'CpxNumTns'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, size(D,3), 'int');
    rD = real(D);    iD = imag(D);
    T = [rD(:)'; iD(:)'];
    fwrite(fid, T, 'double');
   case 'NumVec'
    fwrite(fid, numel(D), 'int');
    for k=1:numel(D)
      serialize(fid, D{k}, type{2});
    end
   case 'NumMat'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    for k=1:numel(D)
      serialize(fid, D{k}, type{2});
    end
   case 'NumTns'
    fwrite(fid, size(D,1), 'int');
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, size(D,3), 'int');
    for k=1:numel(D)
      serialize(fid, D{k}, type{2});
    end
   case 'vectorPoint3'
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, D, 'double');
   case 'vectorIndex3'
    fwrite(fid, size(D,2), 'int');
    fwrite(fid, D, 'int');
   otherwise
    error('wrong');
  end
  
  
  
function setup()
% Adds paths of ButterflyLab to Matlab.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'setup');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'denseSolver']));
addpath(genpath([file_path 'external']));
addpath(genpath([file_path 'matvec']));

end

fprintf('Adding CARLIN directory to PATH\n');
addpath(genpath(pwd));

fprintf('Building MEX\n');
cd('@CARLIN_def');
mex cas9_align_mex.cpp cas9_align.cpp
cd('..');

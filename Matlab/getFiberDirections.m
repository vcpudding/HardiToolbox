function [dirs, weights] = getFiberDirections (x, y, nFibers)
addpath constructMatrixOfIntegrals/
addpath cvx
addpath cvx/structures
addpath cvx/lib
addpath cvx/functions
addpath cvx/commands
addpath cvx/builtins
addpath tensor_toolbox_2.4/
addpath tensor_toolbox_2.4/algorithms/
addpath ~/Libraries/NIFTI_Matlab/

nii = load_nii('dwi-b1500.nii');
grads = dlmread('diffusion_directions.txt');
gradNorms = sum(grads.^2, 2);
GradientOrientations = grads(gradNorms>0, :);

order = 4;
delta = 200;

dwSignal = nii.img(x,y,1,:);
dwSignal = cast(reshape(dwSignal, [numel(dwSignal),1]), 'double');
S = dwSignal(gradNorms>0);
S0 = mean(dwSignal(gradNorms==0));

W=estimate_ODF_with_SDP(S'/S0,GradientOrientations,order,delta);
P=cp_als(tensor(W),nFibers,'init','nvecs','tol',1e-8);

weights = P.lambda/sum(P.lambda);
dirs = P.U{1};

logFile = fopen('log.txt', 'a');
fprintf(logFile, '%d\t%d\t%d\t', x, y, nFibers);
fprintf(logFile, '%f\t', dirs(:));
fprintf(logFile, '\n');
fclose(logFile);

end
% necessary path
addpath(genpath('./Surfacelet'));

% parameters (see mrics5_PLS_TV_init_515_objhalf)
alpha = 2;
alphal = .01;
beta = 1;
betal = .1;
mu = 0.5;
mul = .1;
gamma = .01;

% Can adjust max iteration number and tolerance
MaxIter = 10; 
tol = 1e-4;

% get original P(x,r), choose one slice for experiment, then FFT to S(k,q)
load('Pr_test_12.mat');
slice = 7;          % choose your own slice
Pr = squeeze(Pr(:,:,slice,:,:,:));
Pr = single(500*squeeze(Pr));
N = sqrt(length(Pr(:)));
Skq = single(fftshift(fftn(fftshift(Pr)))/N);

st = tic;
% Choose sample rate and corresponding mask, then apply mask to kq space data
perc = 50;
fname = sprintf('./Syn_msks_12/R_poly_poly_perc%d.mat', perc);
load(fname);
R = logical(R);        
R = fftshift(R);
S = R.*ifftshift(Skq);
mask = R;
            
% zero-filled reconstruction as initialization
recovered_0 = ifftn(S)*N;
Pr_recon0 = real(ifftshift(recovered_0));
Pr_recon0(Pr_recon0<0) = 0;
clear recovered_0

% % run algorithm
tic;
[recovered1, ~, SNR] = mrics5_PLS_TV_init_515_objhalf(Pr_recon0, Pr, R, S, alpha, beta, gamma, mu, alphal, betal, mul, MaxIter, tol, mask);
toc;

% save result           
Pr_recon1 = ifftshift(recovered1);            
Pr_recon1(Pr_recon1<0) = 0;   

% report running time
et = toc(st);
fprintf('RECORD: Running time for sample rate%d is %f', perc, et);

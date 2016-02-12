clear p;
addpath ../
mu = 0; % Parameters for Gaussian
sigma = 1;
% creating p(y|x)
% range for x and y (input and output)
bound_inf = -20; %must be negative (to be chosen)
bound_sup = 20;  %must be positive (to be chosen)
step = 0.1;        % step in quantization (to be chosen)

% E[|X|^rho] < P, rho = 2
interv_P = 1:30;  %choose Power values
rho = 2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interv = bound_inf:step:bound_sup; % support
x = interv;
% p(y|x) = N(y-x), in other words, a shift
% so we can use a toeplitz function to generate p(y|x)
l = interv+abs(bound_inf);
lin = normpdf(x,mu,sigma);  %~ N
p = toeplitz(lin);

optimal = zeros(numel(interv_P),1);
array = zeros(numel(interv_P),1);
for P=interv_P;
   array(P)   = BlahutArimotoConstraint(p,x,P,rho);
   optimal(P) = 1/2*log2(1+P); 
end
figure,plot(array),hold on, plot(optimal,'*')
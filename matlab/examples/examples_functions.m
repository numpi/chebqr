%
% Per compilare il mex file backward_error_pol_mp.c, che fa il conto
% dell'errore all'indietro in multiprecisione, è necessario installare i
% pacchetti libmpfr-dev, libmpc-dev, e poi:
%
% mex backward_error_pol_mp.c -lmpc -lmpfr -lgmp

addpath ..

% Random polynomial, just to test
randpoly100 = chebfun([ randn(100, 1) ; 1], 'coeffs');
randpoly200 = chebfun([ randn(200, 1) ; 1], 'coeffs');
randpoly500 = chebfun([ randn(500, 1) ; 1], 'coeffs');
randpoly1000 = chebfun([ randn(1000, 1) ; 1], 'coeffs');

% Functions to test
fncts = { ...
    % Name of the function -- function handle -- multiprecision required
    { 'p_{100}(x)', randpoly100, true }, ...
    { 'p_{200}(x)', randpoly200, true }, ...
    { 'p_{500}(x)', randpoly500, true }, ...
    { 'p_{1000}(x)', randpoly1000, true }, ...
    { '\log(1 + x + 10^{-3})', @(x) log(1 + 1e-3 + x), true }, ...
    { '\sqrt{x+1.01} - \sin(100x)', @(x) sqrt(1.01 + x) - sin(100*x), true }, ...
    { 'e^x \sin(800x)', @(x) exp(x) .* sin(800 * x), true } , ...
    { '\sin(\frac{1}{x^2 + 10^{-2}})', @(x) sin(1 ./ (x.^2 + 1e-2)), true }, ...
    { 'J_0(20x)', @(x) besselj(0, 20*x), true }, ...
    { '(e^{x^2 - \frac 12} - 1) / (10^{-2} + x^2)', @(x) (exp(x.^2 - .5) - 1) ./ (1e-2 + x.^2), true } , ...
    { '(e^{x^2 - \frac 12} - 1) / (10^{-4} + x^2)', @(x) (exp(x.^2 - .5) - 1) ./ (1e-4 + x.^2), true } , ...
};

fprintf('$f(x)$ & Degree & $\\fastqr$ & \\texttt{eig} & \\texttt{eig\\_nb} \\\\\n \\hline \n')

for j = 1 : length(fncts)
    f = fncts{j}{2};
    fname = fncts{j}{1};
    
    mp_required = fncts{j}{3};
    [be_fastqr, be_eig, be_eignb, deg, p] = test_be_f(f, mp_required);    
    
    fprintf('$%s$ & $%d$ & $%s$ &  $%s$ & $%s$ \\\\ \n', ...
        fname, deg, format_number(be_fastqr), format_number(be_eig), format_number(be_eignb));
end
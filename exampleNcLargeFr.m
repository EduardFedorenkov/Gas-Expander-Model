%% Example. Large Frequency limit
% This example shows the n_cold(r) vs. analytical asymptotic
c = 3 * 10^10;
eVtoErg = 1.6e-12;
mp = 938.27 * 10^6;
mg = mp;
a = 50;
r = linspace(0,a,100);
n0 = 5e+12;
T0 = 0.5;
Twall = 0.026;
np = ones(size(r)) * 1e14;
Tp = ones(size(r)) * 100;
diffCross = 1e-16;
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Nv = 300;

%% Exact solution
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv);

%% Asymptotic solution
eps = 3 / sqrt(2);
alpha = mp / (mg + mp);
VTp = sqrt(2 * Tp / mp) * c;
VT0 = sqrt(2 * T0 / mg) * c;
Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np ./ VTp.^2;
Ueps = eps * VT0 / (2 * alpha) ./ VTp;
nu = 16 * pi * alpha^3 * Q .* VTp.^3 .* exp(-Ueps.^2);
gamma = nu * a / VT0;
x = r/a;
nColdAsym = (n0 / 2) * exp(-(gamma .* (1 - x)).^(2/3)*(2^(1/3) + 2^(-2/3)));

%% Plot 1
figure(1)
plot(x,nCold,'LineWidth',2);
hold on
plot(x,nColdAsym,'LineWidth',2);
lambda = (1 - 5 * VT0 / (nu(1) * a)) * [1, 1];
yLambda = [n0/2, 1];
plot(lambda,yLambda,'g--','LineWidth',2);
legend('Exact solution', 'Asymptotic solution', 'Characteristic length');
lgd = legend;
lgd.FontSize = 18;

%% Plot 2
figure(2)
semilogy(x,nCold,'LineWidth',2);
hold on
semilogy(x,nColdAsym,'LineWidth',2);
legend('Exact solution', 'Asymptotic solution');
lgd = legend;
lgd.FontSize = 18;

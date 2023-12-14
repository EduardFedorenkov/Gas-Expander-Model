c = 3 * 10^10;
np = 1e14;
sigma = 1e-16;
Rp = 50;
T0 = 1;
Tp = 100;
mp = 938.27 * 10^6;
mg = mp;
alpha = mp/(mp+mg);
VT0 = sqrt(2*T0 / mg) * c;
VTp = sqrt(2*Tp / mp) * c;

A = sigma * np * Rp;
a = VTp * 2 * alpha / VT0;

x = linspace(2*alpha,a,3000);
F = A / a * exp(-x.^2 / a^2) ./ x.^2 - exp(-x.^2);
plot(x,F, LineWidth=2);
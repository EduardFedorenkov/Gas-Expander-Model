c = 3 * 10^10;
np = 5e11;
sigma = 1e-16;
Rp = 50;
T0 = 0.5;
Tp = 100;
mp = 938.27 * 10^6;
mg = 2 * mp;
alpha = mp/(mp+mg);
VT0 = sqrt(2*T0 / mg) * c;
VTp = sqrt(2*Tp / mp) * c;

A = sigma * np * Rp;
a = VTp * 2 * alpha / VT0;

x = linspace(0.1,3.5,3000);
y = [-0.9, 0.4];
x_3Sigma = [3/sqrt(2), 3/sqrt(2)];
F1 = A/5 / a * exp(-x.^2 / a^2) ./ x.^2 - exp(-x.^2);
F2 = 0.8 / a * exp(-x.^2 / a^2) ./ x.^2 - exp(-x.^2);
F3 = 3 / a * exp(-x.^2 / a^2) ./ x.^2 - exp(-x.^2);
F4 = 8 / a * exp(-x.^2 / a^2) ./ x.^2 - exp(-x.^2);
plot(x,F1, LineWidth=2);
hold on
plot(x,F2, LineWidth=2);
plot(x,F3, LineWidth=2);
plot(x,F4, LineWidth=2);
plot(x(2171),F1(2171),'r*');
plot(x(1663),F2(1663),'b*');
plot(x(1103),F3(1103),'g*');
plot(x_3Sigma,y,'green');
ylim([-0.9 0.4])
grid on
legend({'A = 0.1','A = 0.8', 'A = 3', 'A = 8', 'root for A = 0.1', ...
    'root for A = 0.8', 'root for A = 3', 'Maxwell distribution 3\sigma'}, ...
    'Location', 'southeast');
xlabel('v/v_{T_0}');
ylabel('\psi(v)');

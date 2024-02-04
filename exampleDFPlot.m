%% Example Title
c = 3 * 10^10;
eVtoErg = 1.6e-12;
mp = 938.27 * 10^6;
mg = 2 * mp;
a = 50;
b = 250;
r = linspace(0,a,35);
R = linspace(a,b,35);
n0 = 5e+12;
T0 = 0.5;
Twall = 0.026;
np = ones(size(r)) * 1e11;
Tp = ones(size(r)) * 100;
diffCross = 1e-16;

Nv = 300;

%% Cold part
% Description of first code block
figure(1);
[DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv, true);

%% Hot part
% Description of second code block
figure(2);
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv);
[DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, Nv, true);

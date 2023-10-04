%Константы задачи
mp = 938.27 * 10^6;
mg = 2 * mp;
a = 50;
r = linspace(0,a,25);
n0 = 6*10^(11);
T0 = 1;
%np = (1 - (r / a).^2) * 1e14;
np = ones(size(r)) * 1e14;
Tp = ones(size(r)) * 2;
diffCross = 1e-16;
%den = ComputeDensityInPlasma(np, Tp, mp, n0, T0, mg, diffCross, r, a);
%plot(r,den); %график концентрации
%DF = PlotDF_R_Phi(np, Tp, mp, n0, T0, mg, diffCross, r, a, 200, 0);
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, 200);
DF = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, 200);
%plot(r,nCold);
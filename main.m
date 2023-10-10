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
%DF = PlotDF_R_Phi(np, Tp, mp, n0, T0, mg, diffCross, r, a, 200, 0);
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, 300);
[DFhot, gridStep] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, 300);
nHot = zeros(1,length(r));
for i = 1:length(r)
    nHot(i) = sum(DFhot(:,:,:,i), "all");
end
nHot = nHot * (gridStep)^3;
figure(26);
plot(r,nCold);
hold on
plot(r,nHot);
plot(r,nCold + nHot);
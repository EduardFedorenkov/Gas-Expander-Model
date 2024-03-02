function [da, IntOnePrr, IntTwoPrr] = AnaliticComputDensePlas(mg, mp, T0, np, Tp, diffCross, a, Nv)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np / VTp^2;
    Ueps = eps * VT0 / (2 * alpha) / VTp;
    nu = 16 * pi * alpha^3 * Q * VTp^3 * exp(-Ueps^2);
    da = 5 * VT0 / nu;
    
    if (a * 0.3 < da)
        disp('Exeption: da seams too large!!!');
    end
    
    VelGridBoundary = eps * VT0 / VTp;                       % V < Veps domen is no need for this solution
    boundaryFactor = zeros(Nv,Nv,Nv);
    
    UGrid = linspace(-1.5 * eps,eps * 1.5,Nv);
    gridStep = UGrid(2) - UGrid(1);
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    UOrtSqrt = sqrt(UOrtSqr);

    N_da = 50;
    factorRange = linspace(0.01,0.3,N_da);           % da / (2 * a)

    IntOnePrr = zeros(size(factorRange));
    IntTwoPrr = zeros(size(factorRange));

    U = zeros(Nv, Nv, Nv);

    for i = 1:Nv
        U(:,:,i) = sqrt(UOrtSqr(:,:) + UGrid(i)^2);
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary^2;
    end

    for j = 1:N_da
        PrrTypeOne = zeros(Nv, Nv, Nv);
        PrrTypeTwo = zeros(Nv, Nv, Nv);

        for i = 1:Nv
            % Prr :
            PrrTypeOne(:,:,i) = exp(-U(:,:,i).^2) ./ (U(:,:,i) .* UOrtSqr) .* Ur.^3 .* ((abs(Uphi) ./ UOrtSqrt) >= (1 - 2 * factorRange(j))) .* (Ur > 0);
            PrrTypeTwo(:,:,i) = exp(-U(:,:,i).^2) ./ U(:,:,i) .* UOrtSqrt .* ((abs(Uphi) ./ UOrtSqrt) < (1 - 2 * factorRange(j))) .* (Ur > 0);

            % qr :
            % PrrTypeOne(:,:,i) = exp(-U(:,:,i).^2) ./ UOrtSqr .* U(:,:,i) .* Ur.^2 .* ((abs(Uphi) ./ UOrtSqrt) >= (1 - 2 * factorRange(j))) .* (Ur > 0);
            % PrrTypeTwo(:,:,i) = exp(-U(:,:,i).^2) ./ Ur .* U(:,:,i) .* UOrtSqrt .* ((abs(Uphi) ./ UOrtSqrt) < (1 - 2 * factorRange(j))) .* (Ur > 0);

            % jr :
            % PrrTypeOne(:,:,i) = exp(-U(:,:,i).^2) ./ (UOrtSqr .* U(:,:,i)) .* Ur.^2 .* ((abs(Uphi) ./ UOrtSqrt) >= (1 - 2 * factorRange(j))) .* (Ur > 0);
            % PrrTypeTwo(:,:,i) = exp(-U(:,:,i).^2) ./ (Ur .* U(:,:,i)) .* UOrtSqrt .* ((abs(Uphi) ./ UOrtSqrt) < (1 - 2 * factorRange(j))) .* (Ur > 0);
        end

        PrrTypeOne = PrrTypeOne .* boundaryFactor;
        PrrTypeTwo = PrrTypeTwo .* boundaryFactor;
    
        IntOnePrr(j) = sum(PrrTypeOne, "all") * gridStep^3;
        IntTwoPrr(j) = sum(PrrTypeTwo, "all") * gridStep^3;
    end

    plot(factorRange, IntOnePrr + IntTwoPrr .* factorRange);
    hold on

    p1 = polyfit(factorRange, IntOnePrr + IntTwoPrr .* factorRange, 1)
    f1 = polyval(p1,factorRange);
    plot(factorRange,f1,'g--');
end


function [DF, gridStep, Vr, Vsqr] = PlotDFCold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv, isPlot)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np ./ VTp.^2;
    Ueps = eps * VT0 / (2 * alpha) ./ VTp;
    nu = 16 * pi * alpha^3 * Q .* VTp.^3 .* exp(-Ueps.^2);
    % Make U grid, Ur and Uphi variables for plot 
    UGrid = linspace(-eps,eps,Nv);
    gridStep = (UGrid(2) - UGrid(1)) * VT0;
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Vr = Ur * VT0;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    Vsqr = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        Vsqr(:,:,i) = UOrtSqr(:,:) + UGrid(i)^2;
    end
    Vsqr = Vsqr * VT0^2;

    VelGridBoundary = eps * VT0;
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) <= VelGridBoundary^2;
    end

    f0 = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        f0(:,:,i) = n0 / (sqrt(pi) * VT0)^3 * exp( - UOrtSqr - Uz(i)^2);
    end
    DF = zeros(Nv, Nv, Nv, length(r));
    for j = 1:length(r)
        for i = 1:Nv
            DF(:,:,i,j) = f0(:,:,i) .* exp(-nu(j) * a / VT0 * (r(j) / a * Ur ./ UOrtSqr  + sqrt(UOrtSqr - (r(j) / a)^2 * Uphi.^2) ./ UOrtSqr));
        end
        DF(:,:,:,j) = DF(:,:,:,j) .* boundaryFactor;
    end
    
    if (isPlot)
        tiledlayout(3,3);

        ax1 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),1));
        shading flat;
        shading interp;
        title(ax1,'f_g(r = 0)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\phi/v_{Tg}');
    
        ax2 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),2));
        shading flat;
        shading interp;
        title(ax2,'f_g(r = 0.125a)'); 
        xlabel(ax2,'v_r/v_{Tg}'); 
        ylabel(ax2,'v_\phi/v_{Tg}');
    
        ax3 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),3));
        shading flat;
        shading interp;
        title(ax3,'f_g(r = 0.25a)'); 
        xlabel(ax3,'v_r/v_{Tg}'); 
        ylabel(ax3,'v_\phi/v_{Tg}');
    
        ax4 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),4));
        shading flat;
        shading interp;
        title(ax4,'f_g(r = 0.375a)'); 
        xlabel(ax4,'v_r/v_{Tg}'); 
        ylabel(ax4,'v_\phi/v_{Tg}');
    
        ax5 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),5));
        shading flat;
        shading interp;
        title(ax5,'f_g(r = 0.5a)'); 
        xlabel(ax5,'v_r/v_{Tg}'); 
        ylabel(ax5,'v_\phi/v_{Tg}');
    
        ax6 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),6));
        shading flat;
        shading interp;
        title(ax6,'f_g(r = 0.625a)'); 
        xlabel(ax6,'v_r/v_{Tg}'); 
        ylabel(ax6,'v_\phi/v_{Tg}');
    
        ax7 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),7));
        shading flat;
        shading interp;
        title(ax7,'f_g(r = 0.75a)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\phi/v_{Tg}');
    
        ax8 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),8));
        shading flat;
        shading interp;
        title(ax8,'f_g(r = 0.875a)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\phi/v_{Tg}');
    
        ax9 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),9));
        shading flat;
        shading interp;
        title(ax9,'f_g(r = a)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\phi/v_{Tg}');
    end
end

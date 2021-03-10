function [V0,DomS,Domc,Som] = PublicPerturbativeHigherHarmonicPotential(Radius,dr,NHarmonics,Vcoeff,omegaMax,S1,omega,LinRef)
%{

FUNCTION DESCRIPTION
    
    Generates the matrix operator which acts on the linear state
    vectors S = [S3 S5 ... ] and C = [C3 C5 ... ], and the source term.
    
=====================================================================

INPUT DESCRIPTION

    Radius = Physical size of the grid. If Radius is chosen too large, then
    there may not be enough precision to satisfy the boundary conditions.

    dr = Physical size of the lattice spacing. 

    NHarmonics = Number of harmonics to calculate.

    Vcoeff = The potential coefficients, whose sum must be less than or
    equal to omegaMax^2.

    S1 = Oscillon fundamental mode.

    omega = Fundamental frequency of the oscillon in units of the mass.
    
    LinRef = Number of additional grid-points per grid-point of the S1
    shooting code.
    
=====================================================================

%}
    
    if sum(Vcoeff) ~= omegaMax^2
        Vcoeff = [Vcoeff omegaMax^2 - sum(Vcoeff)];
    end
    
    GridSize = length(S1);
    BigGrid = LinRef * GridSize;
    r = linspace(0,Radius,GridSize)';
    dr = dr/LinRef;
    rHR = linspace(0,Radius,BigGrid)';
    S1 = interp1(r,S1,rHR);
    
    V0 = zeros(BigGrid * NHarmonics,1);
    DomSContainer = zeros(BigGrid * NHarmonics, 2 * NHarmonics - 1);
    DomcContainer = zeros(BigGrid * NHarmonics, 2 * NHarmonics - 1);
    
    NabC = vertcat(ones(BigGrid - 1,1), dr/2);
    OmegaMultiplier = vertcat(ones(BigGrid - 1,1),0);
    
    for harmonicN = 1 : NHarmonics
        for harmonicL = 1 : NHarmonics
            N = 2 * harmonicN + 1;
            L = 2 * harmonicL + 1;
            VNLS = zeros(BigGrid,1);
            VNLc = zeros(BigGrid,1);
            if N == L
                VNLS = VNLS - (-2/(dr^2)) * NabC - ((N * omega)^2) * OmegaMultiplier;
                VNLc = VNLc - (-2/(dr^2)) * NabC - ((N * omega)^2) * OmegaMultiplier;
            end
            for mm = 1 : length(Vcoeff)
                VNLS = VNLS +  (Vcoeff(mm)/(omegaMax^2)) * (besselj(N - L,mm * S1/omegaMax) - besselj(N + L,mm * S1/omegaMax)) .* OmegaMultiplier;
                VNLc = VNLc +  (Vcoeff(mm)/(omegaMax^2)) * (besselj(N - L,mm * S1/omegaMax) + besselj(N + L,mm * S1/omegaMax)) .* OmegaMultiplier;
            end
            row = harmonicN;
            column = NHarmonics + harmonicN - harmonicL;
            DomSContainer((row - 1) * BigGrid + 1:row * BigGrid,column) = VNLS;
            DomcContainer((row - 1) * BigGrid + 1:row * BigGrid,column) = VNLc;
        end
    end
    
    NabL = -1/(dr^2) * ones(BigGrid * NHarmonics,1);
    NabR = NabL;
    
    for harmonic = 1 : NHarmonics - 1
        NabL(BigGrid * harmonic) = 0;
        NabR(BigGrid * harmonic + 1) = 0;
        NabR(BigGrid * harmonic) = -1/dr;
    end
    NabR(BigGrid * NHarmonics) = -1/dr;
    
    DomSContainer = [DomSContainer(:,1 : NHarmonics - 1) NabL DomSContainer(:,NHarmonics) NabR DomSContainer(:,NHarmonics + 1 : 2 * NHarmonics - 1)];
    DomcContainer = [DomcContainer(:,1 : NHarmonics - 1) NabL DomcContainer(:,NHarmonics) NabR DomcContainer(:,NHarmonics + 1 : 2 * NHarmonics - 1)];
    
    LeftIdx = -(NHarmonics - 1) * BigGrid:BigGrid : - BigGrid;
    RightIdx = BigGrid : BigGrid : (NHarmonics - 1) * BigGrid;
    SparseDiagIndex = [LeftIdx (-1) 0 (1) RightIdx];
    
    DomS = spdiags(DomSContainer, SparseDiagIndex, BigGrid * NHarmonics, BigGrid * NHarmonics);
    Domc = spdiags(DomcContainer, SparseDiagIndex, BigGrid * NHarmonics, BigGrid * NHarmonics);
    
    
    
    SomMultiplier = vertcat(zeros(BigGrid - 2,1),0,1);
    SomContainer = zeros(BigGrid * NHarmonics,1);
    for harmonicN = 1 : NHarmonics
        N = 2 * harmonicN + 1;
        SomContainer((harmonicN - 1) * BigGrid + 1:harmonicN * BigGrid) = sqrt((N * omega)^2 - 1) * SomMultiplier;
    end
    
    Som = spdiags(SomContainer, 0, BigGrid * NHarmonics, BigGrid * NHarmonics);
    
    
    for harmonic = 1 : NHarmonics
        
        VN = zeros(BigGrid,1);
        N = 2 * harmonic + 1;
        
        for mm = 1 : length(Vcoeff)
            VN = VN + 2 * (Vcoeff(mm)/(mm * omegaMax)) * besselj(N,mm * S1/omegaMax) .* rHR;
        end
        V0((harmonic - 1) * BigGrid + 1 : harmonic * BigGrid) = VN;
    end
end

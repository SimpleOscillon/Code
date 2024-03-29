function [S, C, r] = PublicPerturbativeOscillon(Radius,dr,Vcoeff,thetaMax,NHarmonics,omega,LinRef,S10,NIterations)
%{

FUNCTION DESCRIPTION

    Solves the linear set of equations describing the evolution of the
    higher harmonics in the background of the fundamental mode S1. 

=====================================================================


INPUT DESCRIPTION

    Radius = Physical size of the S1 grid. If Radius is chosen too large, then
    there may not be enough precision to satisfy the boundary conditions.

    dr = Physical size of the lattice spacing. 

    thetaMax = periodicity of the potential.

    Vcoeff = The potential coefficients, whose sum must be less than or
    equal to omegaMax^2.

    NHarmonics = Number of harmonics to calculate.

    omega = Physical frequency of the oscillon.
    
    LinRef = Number of additional grid-points added to the linear harmonic equation
    per grid point in the S1 shooting code.

    S10 = One tenth of the range of central values to scan over in the shooting code.
=====================================================================

%}

    if sum(Vcoeff) ~= thetaMax^2
        Vcoeff = [Vcoeff thetaMax^2 - sum(Vcoeff)];
    end
    
    HigherHarmonicS = [];
    
    for iteration = 1 : NIterations
    
        S1 = PublicFundamentalModeShooting(Radius,dr,Vcoeff,thetaMax,omega,S10,HigherHarmonicS,NHarmonics);

        [V0,DomS,Domc,Som] = PublicPerturbativeHigherHarmonicPotential(Radius,dr,NHarmonics,Vcoeff,thetaMax,S1,omega,LinRef);

        BigMatrix = [DomS -Som;Som Domc];
        BigVector = [-V0;zeros(size(V0))];
               
        SC = reshape(BigMatrix\BigVector,length(S1) * LinRef,2 * NHarmonics);
        Stemp = [S1 SC(1:LinRef:end,1 : NHarmonics)];
        Ctemp = [zeros(size(S1)) SC(1:LinRef:end,NHarmonics + 1 : 2 * NHarmonics)];
        SReCompute = [zeros(size(S1)) SC(1:LinRef:end,1 : NHarmonics) zeros(size(S1))];
        
        r = linspace(0,Radius,length(S1))';
        HigherHarmonicS = SReCompute./(r + dr * (1/LinRef));
    end
    
    GridSize = length(S1) - 1;
    r = linspace(0,Radius,GridSize)' + dr * (1/LinRef);
    S = Stemp(1:end-1,:);
    C = Ctemp(1:end-1,:);
    
end

function [S, C, r] = PublicPerturbativeOscillon(Radius,RadiusRad,dr,Vcoeff,omegaMax,NHarmonics,omega,LinRef,S10)
%{

FUNCTION DESCRIPTION

    Solves the linear set of equations describing the evolution of the
    higher harmonics in the background of the fundamental mode S1. 

=====================================================================


INPUT DESCRIPTION

    Radius = Physical size of the S1 grid. If Radius is chosen too large, then
    there may not be enough precision to satisfy the boundary conditions.

    RadiusRad = Physical size of the grid on which the radiation is
    calculated. Generally, we must have RadiusRad ≥ Radius. Outside of
    Radius, S1 is taken to be identically zero.

    dr = Physical size of the lattice spacing. 

    Vcoeff = The potential coefficients, whose sum must be less than or
    equal to 1.

    NHarmonics = Number of harmonics to calculate.

    omega = Physical frequency of the oscillon.


=====================================================================

%}

    if sum(Vcoeff) ~= omegaMax^2
        Vcoeff = [Vcoeff omegaMax^2 - sum(Vcoeff)];
    end
    
    RadiusRad = RadiusRad + dr;
    
    NRad = round(RadiusRad/dr);
    
    S1 = PublicFundamentalModeShooting(Radius,dr,Vcoeff,omegaMax,omega,S10);
    S1 = vertcat(S1,zeros(NRad - length(S1),1));
    
    [V0,DomS,Domc,Som] = PublicPerturbativeHigherHarmonicPotential(RadiusRad,dr,NHarmonics,Vcoeff,omegaMax,S1,omega,LinRef);
    
    BigMatrix = [DomS -Som;Som Domc];
    BigVector = [-V0;zeros(size(V0))];
    
    SC = reshape(BigMatrix\BigVector,length(S1) * LinRef,2 * NHarmonics);
    Stemp = [S1 SC(1:LinRef:end,1 : NHarmonics)];
    Ctemp = [zeros(size(S1)) SC(1:LinRef:end,NHarmonics + 1 : 2 * NHarmonics)];
    S = Stemp(1:end-1,:);
    C = Ctemp(1:end-1,:);
    GridSize = length(S1) - 1;
    r = linspace(0,Radius,GridSize)';
    
end
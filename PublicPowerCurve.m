function [PowerVsOmegaList,EnergyVsOmegaList,Lifetime,PowerListHarmonics,SList,CList,r] = PublicPowerCurve(Radius,dr,Vcoeff,omegaMax,NHarmonics,OmegaList,LinRef,S10)
%{

FUNCTION DESCRIPTION

    Treating only the fundamental mode as non-perturbative, calculates the
    power radiated by the higher harmonics if the configuration is
    linearly stable. Stores the output as PowerVsOmega = [OmegaList
    PowerList].

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

    OmegaList = List of the physical frequencies of the oscillon at which
    to compute the output power.


=====================================================================

%}
    RadiusRad = Radius; %The radius out to which the perturbative harmonics are computed.
                        %This can be changed to a larger value.
                        
    PowerListHarmonics = zeros(length(OmegaList),NHarmonics);
    EnergyList = zeros(size(OmegaList));
    MassList = zeros(size(OmegaList));
    SList = [];
    CList = [];
    
    if sum(Vcoeff) ~= omegaMax^2
        Vcoeff = [Vcoeff omegaMax^2 - sum(Vcoeff)];
    end
    
    index = 0;
    for omega = OmegaList
        disp(omega)
        index = index + 1;
        r = linspace(0,Radius,(Radius/dr) + 1)';
        [S, C, ~] = PublicPerturbativeOscillon(Radius,RadiusRad,dr,Vcoeff,omegaMax,NHarmonics,omega,LinRef,S10);
        
        SList = [SList S];
        CList = [CList C];
        
        SVal = S(end,:);
        CVal = C(end,:);
        FrequencyCoefficient = 1 : 2 : 2 * NHarmonics + 1;
        
        PowerListHarmonics(index,:) = -2 * pi * FrequencyCoefficient(2:end) .* sqrt(omega^2 * FrequencyCoefficient(2:end).^2 - 1) .* (SVal(2:end).^2 + CVal(2:end).^2);
        
        S1Derivative = (S(2 : end,1)-S(1 : end-1,1))/dr;
        S1Val = S(1:end-1,1);
        S1Potential = 0;
        for mm = 1 : length(Vcoeff)
            S1Potential = S1Potential + (Vcoeff(mm)/mm^2) * (1 - besselj(0,mm * S1Val/omegaMax));
        end
        r = r(1 : length(S1Val));
        EnergyList(index) = sum(4 * pi * dr * r.^2 .* ((S1Derivative/2).^2 + omega^2 * (S1Val/2).^2 + S1Potential));
        MassList(index) = sum(4 * pi * dr * r.^2 .*((S1Val/2).^2));
        
    end
    
   
    
    PowerList = sum(PowerListHarmonics,2);
    dEdOmega = (EnergyList(1:end-1) - EnergyList(2:end));
    ViableFrequencies = dEdOmega > 0;

    Lifetime = sum(- ViableFrequencies .*dEdOmega./ PowerList(2:end)');
    PowerVsOmegaList = [OmegaList' PowerList];
    EnergyVsOmegaList = [OmegaList' EnergyList'];
    
end

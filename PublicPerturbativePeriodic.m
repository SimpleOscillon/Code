
format long

                %The Fourier coefficients of the potential, normalized to
                %sum to omegaMax^2. If the total of Vcoeff does not sum to
                %omegaMax^2, then it an additional term will automatically
                %be added to satisfy the mass constraint.
Vcoeff = [];


                %Fundamental periodicity of the periodic potential in units
                %of f.
omegaMax = 1;   


Radius = 15;    %Radius out to which the fundamental bound harmonic is computed.
dr = 0.01;
LinRef = 10;    %Number of additional grid points in the radiation computation
                %per grid point in the fundamental mode. Higher LinRef improves
                %the resolution of the boundary condition at the origin.
                
S10 = 5;        %Shooting range: if result does not converge make this number
                %larger or smaller.

                
               
                

NHarmonics = 2; %Number of perturbative harmonics to compute
OmegaList = 0.80:0.01:0.94;
[PowerVsOmegaList,EnergyVsOmegaList,Lifetime,PowerInHarmonics,SList,CList,r]...
    = PublicPowerCurve(Radius,dr,Vcoeff,omegaMax,NHarmonics,OmegaList,LinRef,S10);
                                            %these last four arguments can
                                            %be used to extract the 
                                            %1) Power stored in each
                                            %harmonic,
                                            %2) The PQB profile, multiplied
                                            %by the radius r,
                                            %3) The OD profile, also
                                            %multiplied by the radius r,
                                            %4) The corresponding list of
                                            %radius coordinates, r.
dE = EnergyVsOmegaList(2:end,2) - EnergyVsOmegaList(1:end-1,2) ;

log10(Lifetime)
figure(1)
hold on
plot(PowerVsOmegaList(2:end,1),log10(-PowerVsOmegaList(2:end,2).* (dE < 0)))
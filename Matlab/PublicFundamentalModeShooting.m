function S1 = PublicFundamentalModeShooting(Radius,dr,Vcoeff,thetaMax,omega,S10,HigherHarmonicS,NHarmonics)

%{

FUNCTION DESCRIPTION

    Returns the fundamental mode profile of the PQB, neglecting the
    contribution of higher harmonics.

=====================================================================


INPUT DESCRIPTION
    
    Radius = Physical size of the grid. If Radius is chosen too large, then
    there may not be enough precision to satisfy the boundary conditions.

    dr = Physical size of the lattice spacing. 

    Vcoeff = The potential coefficients, whose sum must be less than or
    equal to thetaMax^2.
    
    thetaMax = Fundamental periodicity of the potential.

    omega = The physical frequency of the oscillon.
    
    S10 = One tenth of the range of central field amplitudes to search over 
    in the shooting code. If the code starts outputting nonsense, try 
    adjusting this.

=====================================================================

%}


    nextidx = -1;
    
    %Determine the size of the grid.
    gridSize = round(Radius / dr);
    
    %Enforce that Vcoeff sums to omegaMax^2.
    if sum(Vcoeff) ~= thetaMax^2
        Vcoeff = [Vcoeff thetaMax^2 - sum(Vcoeff)];
    end
    
    if isempty(HigherHarmonicS) == 1
    
        %Number of digits known about initial condition.
        for digit = 0 : 15

            %Set of initial conditions to check.
            init = ones(19,1)* S10 + (-9:1:9)' * 10^(-digit);
            positive = zeros(19,1);

            %Looping over the set of initial conditions.
            for jj = 1 : 19


                S1 = zeros(gridSize + 1,1);

                S1(1) = init(jj);
                S1(2) = init(jj);

                %Evolve from the initial condition
                for ii = 1 : gridSize - 1

                    %Calculate the effective potential
                    Veff = 0;
                    for mm = 1 : length(Vcoeff)
                        Veff = Veff + (Vcoeff(mm)/(mm * thetaMax)) * besselj(1,mm * S1(ii + 1)/thetaMax);
                    end

                    %Finite difference solver
                    S1(ii + 2) = (2 * dr^2) * Veff + (-1 + (2/ii)) * S1(ii) + (2 - (2/ii)  - (dr * omega)^2) * S1(ii + 1);

                end

                positive(jj) = min(S1);
                lastidx = nextidx;
                nextidx = find(positive > 0, 1, 'last');

                if lastidx == nextidx
                    break
                end

            end
        if isempty(nextidx) && digit ~= 0
            break
        elseif isempty(nextidx) && digit == 0
            S10 = S10/10;
        elseif isempty(nextidx) == false
            S10 = round(init(nextidx),digit);
        end

        end
    else
        %Number of digits known about initial condition.
    for digit = 0 : 15
    
        %Set of initial conditions to check.
        init = ones(19,1)* S10 + (-9:1:9)' * 10^(-digit);
        positive = zeros(19,1);

        %Looping over the set of initial conditions.
        for jj = 1 : 19

            
            S1 = zeros(gridSize + 1,1);

            S1(1) = init(jj);
            S1(2) = init(jj);
            
            %Evolve from the initial condition
            for ii = 1 : gridSize - 1
                
                %Calculate the effective potential
                Veff = 0;
                for mm = 1 : length(Vcoeff)
                    for harmonic = 1 : NHarmonics
                        harmonicNumber = 2 * harmonic - 1;
                        if harmonic == 1
                            Veff = Veff + 2 * (Vcoeff(mm)/(mm * thetaMax)) * besselj(1,mm * S1(ii + 1)/thetaMax);
                        else
                            Veff = Veff + (Vcoeff(mm)/(thetaMax^2)) * besselj(harmonicNumber - 1,mm * S1(ii + 1)/thetaMax) ...
                                * (HigherHarmonicS(ii + 1,harmonic) - HigherHarmonicS(ii + 1,harmonic - 1))/(ii * dr);
                        end
                    end
                end
                
                %Finite difference solver
                S1(ii + 2) = (dr^2) * Veff + (-1 + (2/ii)) * S1(ii) + (2 - (2/ii)  - (dr * omega)^2) * S1(ii + 1);
                
            end
   
            positive(jj) = min(S1);
            lastidx = nextidx;
            nextidx = find(positive > 0, 1, 'last');
            
            if lastidx == nextidx
                break
            end
            
        end
    if isempty(nextidx) && digit ~= 0
        break
    elseif isempty(nextidx) && digit == 0
        S10 = S10/10;
    elseif isempty(nextidx) == false
        S10 = round(init(nextidx),digit);
    end
        
    end
end

function S1 = PublicFundamentalModeShooting(Radius,dr,Vcoeff,omegaMax,omega,S10)

%{

FUNCTION DESCRIPTION

    Returns the fundamental mode profile of the oscillon, neglecting the
    contribution of higher harmonics.

=====================================================================


INPUT DESCRIPTION
    
    Radius = Physical size of the grid. If Radius is chosen too large, then
    there may not be enough precision to satisfy the boundary conditions.

    dr = Physical size of the lattice spacing. 

    Vcoeff = The potential coefficients, whose sum must be less than or
    equal to 1.

    omega = Physical frequency of the oscillon.

=====================================================================

%}


    nextidx = -1;
    
    %Determine the size of the grid.
    gridSize = round(Radius / dr);
    
    %Enforce that Vcoeff sums to 1.
    if sum(Vcoeff) ~= omegaMax^2
        Vcoeff = [Vcoeff omegaMax^2 - sum(Vcoeff)];
    end
    
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
                    Veff = Veff + (Vcoeff(mm)/(mm * omegaMax)) * besselj(1,mm * S1(ii + 1)/omegaMax);
                end
                
                %Finite difference solver
                S1(ii + 2) = (2 * dr^2) * Veff + (-1 + (2/ii)) * S1(ii) + (2 - (2/ii)  - (dr * omega)^2) * S1(ii + 1);
                
            end
%             
%             figure(1)
%             hold on
%             plot(S1)
%             pause(0.1)
%             
            positive(jj) = min(S1);
            
            %Find the smallest initial condition that leads to an always
            %positive S1
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
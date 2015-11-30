function [diag, sub, sup, rhs] = Assemble_BeamWarming( u, epsilon, dt, dx )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the Beam and Warming system
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, November 2015
    %%%
    
    F = u.^2 / 2;
    
    N = length(u);
    
    diag_range = 2:N-1;
     sub_range = 3:N-1;
     sup_range = 2:N-2;
    
    diag = ones(length(diag_range),1);
     sub = - (1/4) * (dt/dx) *  u(sub_range-1);
     sup =   (1/4) * (dt/dx) *  u(sup_range+1);
     rhs = u(diag_range) ...
           - (1/2) * (dt/dx) * (F(diag_range+1) - F(diag_range-1)) ...
           + (1/4) * (dt/dx) *  u(diag_range+1).^2 ...
           - (1/4) * (dt/dx) *  u(diag_range-1).^2 ;
    
    % Account for boundaries.
    rhs(1)   = rhs(1)   + (1/4) * (dt/dx) * 10;
    rhs(end) = rhs(end) - (1/4) * (dt/dx) *  0;
    
    % Add the artificial viscosity.
    Deps = zeros(length(diag_range),1);
    for i = 1:length(diag_range)
        ii = i+1;
        
        if ii-2 > 0
            Deps(i) = Deps(i) +     u(ii-2);
        else
            Deps(i) = Deps(i) +     10;
        end
        
        if ii-1 > 0
            Deps(i) = Deps(i) - 4 * u(ii-1);
        else
            Deps(i) = Deps(i) - 4 * 10;
        end
        
            Deps(i) = Deps(i) + 6 * u(ii);
            
        if ii+1 < N+1
            Deps(i) = Deps(i) - 4 * u(ii+1);
        else
            Deps(i) = Deps(i) - 4 *  0;
        end
        
        if ii+2 < N+1
            Deps(i) = Deps(i) +     u(ii+2);
        else
            Deps(i) = Deps(i) +      0;
        end
    end
    Deps = -epsilon * Deps;
    rhs = rhs + Deps;

end
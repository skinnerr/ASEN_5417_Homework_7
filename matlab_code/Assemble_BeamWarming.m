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
    
    % Add dissipation.
    De = zeros(length(diag_range),1);
    for i = 1:length(diag_range)
        ii = i+1;
        
        if ii-2>0 tmp=u(ii-2);      else tmp=10; end
        De(i) = De(i) + tmp;
        
        if ii-1>0 tmp = -4*u(ii-1); else tmp=-4*10; end
        De(i) = De(i) + tmp;
        
        De(i) = De(i) + 6 * u(ii);
            
        if ii+1<N+1 tmp=-4*u(ii+1); else tmp=0; end
        De(i) = De(i) + tmp;
        
        if ii+2<N+1 tmp=u(ii+2);    else tmp=0; end
        De(i) = De(i) + tmp;
    end
    De = -epsilon * De;
    rhs = rhs + De;

end

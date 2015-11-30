function [] = Problem_1()

    %%%%%%
    % Solves the inviscid Burgers equation using the MacCormack explicit method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Solution domain.
    dx = 1;
    x = (0:dx:40)';
    N = length(x);
    t = 0;
    
    % Initialize the solution, indexed by (x,y), and set BCs.
    u = zeros(N,1);
    for i = 1:length(x)
        if x(i) <= 30
            u(i) = 10;
        end
    end
    
    %%%
    % Solve problem numerically.
    %%%
    
    % For each Courant number...
%     for C = [0.75, 1.00, 1.25]
    for C = 0.75

        % Solve problem up to t~8.
        while t(end) < 8
            
            u_n = u(:,end);
            u_hat = nan(N,1);
            u_np1 = nan(N,1);

            % Calculate time step.
            dt = C * dx / max(u(:,end));
            t(end+1) = t(end) + dt;
            
            % MacCormack iterations.
            for i = 2:N
                u_hat(i) = u_n(i) - u_n(i) * (dt/dx) * (u_n(i) - u_n(i-1));
            end
            u_hat(1) = u_hat(2)
            for i = 1:N-1
                u_np1(i) = 0.5 * ( (u_n(i) + u_hat(i)) ...
                                  - u_n(i) * (dt/dx) * (u_hat(i+1) - u_hat(i)));
            end
            u_np1(N) = u_np1(N-1)
            
            % Update solution.
            u(:,end+1) = u_np1;

        end
        
    end

    %%%
    % Process results.
    %%%
    figure();
    plot(t);
    figure();
    surf(x,t,u');
    xlabel('x');
    ylabel('t');
    xlim([0,40]);
    ylim([min(t),max(t)]);
    
    disp('Done.');
    return
    
end














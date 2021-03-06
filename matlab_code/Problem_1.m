function [] = Problem_1()

    %%%%%%
    % Solves the inviscid Burgers equation using the MacCormack explicit method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % For each Courant number...
    for C = [0.75, 1.00, 1.25]
        
        %%%
        % Define variables specific to the boundary-value problem.
        %%%

        % Solution domain.
        dx = 0.2;
        x_min = 0;
        x_max = 60;
        x = (x_min:dx:x_max)';
        N = length(x);
        t = 0;
        dt_history = 0;

        % Initialize the solution, indexed by (x,y), and set BCs.
        u = zeros(N,1);
        for i = 1:length(x)
            if x(i) <= 15
                u(i) = 10;
            end
        end

        %%%
        % Solve problem numerically.
        %%%

        % Solve problem up to t~8.
        while t(end) < 8
            
            u_n = u(:,end);
            F_n = u_n.^2 / 2;
            u_h = nan(N,1);
            u_np1 = nan(N,1);

            % Calculate time step.
            dt = C * dx / max(u_n);
            dt_history(end+1) = dt;
            t(end+1) = t(end) + dt;
            
            % MacCormack iterations.
            
            calc_u_hat = @(uni,Fnip1,Fni)       uni                                  ...
                                                    - (dt/dx) * (Fnip1 - Fni);
            calc_u_np1 = @(uni,uhi,Fhip1,Fhi) ( uni + uhi                            ...
                                                    - (dt/dx) * (Fhip1 - Fhi) ) / 2;
                
            for i = 1:N-1
                u_h(i) = calc_u_hat(u_n(i), F_n(i+1), F_n(i));
            end
            u_h(N) = u_h(N-1);
            F_h = u_h.^2 / 2;

            for i = 2:N
                u_np1(i) = calc_u_np1(u_n(i), u_h(i), F_h(i), F_h(i-1));
            end
            u_np1(1) = u_np1(2);
            
            % Update solution.
            u(:,end+1) = u_np1;

        end

        %%%
        % Process results.
        %%%
        
        % Time evolution of solution for a single Courant number.
        n_plot = 5;
        cmap = winter(n_plot);
        step_numbers = round(linspace(1,length(t),n_plot));
        hf = figure(round(C*10));
        set(hf,'Position',[100,500,900,300]);
        hold on;
        for i = 1:length(step_numbers)
            tmp = sprintf('t = %.2f', t(step_numbers(i)));
            plot(x, u(:,step_numbers(i)), 'DisplayName', tmp, 'Color', cmap(i,:));
        end
        title(sprintf('C = %.2f',C));
        xlabel('x');
        ylabel('u');
        ylim([0,15]);
        xlim([x_min,x_max]);
        hleg = legend('show');
        set(hleg,'Location','southwest');
        
        % Solution at t~8, comparing Courant numbers.
        hf = figure(1);
        set(hf,'Position',[100,500,900,300]);
        hold on;
        tmp = sprintf('C = %.2f, t = %.2f', C, t(end));
        plot(x, u(:,end), 'DisplayName', tmp);
        xlabel('x');
        ylabel('u');
        ylim([0,15]);
        xlim([50,60]);
        
    end
    
    figure(1);
    hleg = legend('show');
    set(hleg,'Location','southwest');
    
    disp('Done.');
    return
    
end














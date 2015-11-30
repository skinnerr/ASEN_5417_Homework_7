function [] = Problem_2()

    %%%%%%
    % Solves the inviscid Burgers equation using the Beam and Warming implicit method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % For each Courant number...
    for C = [0.75, 1.00, 1.25]
%     for C = 1.00
        
        %%%
        % Define variables specific to the boundary-value problem.
        %%%

        % Solution domain.
        dx = 0.2;
        x = (0:dx:40)';
        N = length(x);
        t = 0;
        dt_history = 0;

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
        
        direction_flag = true;

        % Solve problem up to t~8.
        while t(end) < 8
            
            u_n = u(:,end);
            u_np1 = nan(N,1);
            dt = 0.1;
            t(end+1) = t + dt;
            
            % Beam and Warming implicit iteration.
            for i = 2:N-1
                [diag, sub, sup, rhs] = Assemble_BeamWarming(u_n, );
                [sol] = Thomas(diag, sub, sup, rhs);
                u(:,j) = [BC.uw; sol; BC.ue];
            end

            if any(u_np1 > 1e3)
                fprintf('Courant number %.2f diverged.', C);
                break
            end
            
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
        set(hf,'Position',[100,500,1100,350]);
        hold on;
        for i = 1:length(step_numbers)
            tmp = sprintf('t = %.2f', t(step_numbers(i)));
            plot(x, u(:,step_numbers(i)), 'DisplayName', tmp, 'Color', cmap(i,:));
        end
        title(sprintf('C = %.2f',C));
        xlabel('x');
        ylabel('u');
        ylim([0,12.5]);
        xlim([0,40]);
        hleg = legend('show');
        set(hleg,'Location','southwest');
        
        % Solution at t~8, comparing Courant numbers.
        hf = figure(1);
        set(hf,'Position',[100,500,1100,350]);
        hold on;
        tmp = sprintf('C = %.2f, t = %.2f', C, t(end));
        plot(x, u(:,end), 'DisplayName', tmp);
        xlabel('x');
        ylabel('u');
        ylim([0,12.5]);
        xlim([0,40]);
        
        figure();
        surf(x,t,u');
        xlabel('x');
        ylabel('t');
        title(sprintf('C = %.2f',C));
        xlim([0,40]);
        ylim([min(t),max(t)]);
        zlim([0,15]);
        
    end
    
    figure(1);
    hleg = legend('show');
    set(hleg,'Location','southwest');
    
    disp('Done.');
    return
    
end










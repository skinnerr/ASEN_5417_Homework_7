function [] = Problem_2()

    %%%%%%
    % Solves the inviscid Burgers equation using the Beam and Warming implicit method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    clear
    
    % For each Courant number...
    for C = [0.75, 1.00, 1.25]
        epsilon = 0.1;
        
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
            dt = C * dx / max(u_n);
            t(end+1) = t(end) + dt;
            
            % Beam and Warming implicit iteration.
            [diag, sub, sup, rhs] = Assemble_BeamWarming(u_n, epsilon, dt, dx);
            [sol] = Thomas(diag, sub, sup, rhs);
            u_np1 = [10; sol; 0];
            u_np1(1:10) = 10;
            
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
        text(50.1,2.5,sprintf('Epsilon = %.3f', epsilon),'FontSize',14);
        xlabel('x');
        ylabel('u');
        ylim([0,15]);
        xlim([50,60]);
        
    end

    figure(1);
    hleg = legend('show');
    set(hleg,'Location','northeast');
    
    disp('Done.');
    return
    
end














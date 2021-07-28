function Y_vax = gdp_vax(I_grid, R_grid, PAR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GDP_VAX  Compute the continuation value of GDP after arrival of the
    % vaccine. Called by simulations.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   I_grid = initial share of infected individuals.
    %   R_grid = initial share of recovered individuals.
    %   PAR = structure array of parameters.
    %
    % Outputs:
    %   Y_vax = continuation value of GDP after arrival of the vaccine.

    
    %%% initialization
    % setup grids of health state and lockdown policy
    [Ig, Rg] = meshgrid(I_grid, R_grid); % set up grid with I and R. Here the dim of grid is n_R*n_I

    [~, I_prime, R_prime, ~] = sir_model(Ig, Rg, 1, PAR, true); % next period number of each variety of people

    % functional forms
    A = PAR.TFP1;
    S = @(I, R) 1 - I - R*(1 + PAR.RHO4/PAR.RHO3); % share susceptible (S) from rearranging SIR model
    l_star = @(I, R) S(I, R) + PAR.GAMMA*I + R; % effective labor (\ell) incl. lockdown (1 - L)
    x_star = @(I, R) (PAR.ALPHA.*A/PAR.INTEREST_W).^(1/(1 - PAR.ALPHA)).*l_star(I, R); % investment (x) from no-arbitrage condition
    y_star = @(I, R) A.*x_star(I, R).^(PAR.ALPHA).*l_star(I, R).^(1 - PAR.ALPHA); % output (y) from functional form assumption
    
    disc_factor = 1/(1 + PAR.INTEREST_W);
    
    %%% value function iteration
    % set initial value, assume to be SS value
    Y_vax = y_star(Ig, Rg)/(1 - disc_factor); 

    n_I = length(I_grid);
    n_R = length(R_grid);

    yinterp = y_star(Ig(:), Rg(:)) + disc_factor*interp2(I_grid, R_grid, Y_vax, I_prime(:), R_prime(:), 'linear'); 
    yinterp = reshape(yinterp, n_R, n_I); % update

    crit = 1e-8; % critical value for 
    
    % iteration
    while max(max(abs(yinterp - Y_vax))) >= crit
        Y_vax = yinterp;
        yinterp = y_star(Ig(:), Rg(:)) + disc_factor*interp2(I_grid, R_grid, Y_vax, I_prime(:), R_prime(:), 'linear');
        yinterp = reshape(yinterp, n_R, n_I);
    end

    % mute infeasible cases
    for i = 1:1:n_R
        for j = 1:1:n_I
            s = 1 - I_grid(j) - R_grid(i)*(1 + PAR.RHO4/PAR.RHO3);
            if s < 0
               Y_vax(i, j) = -Inf; % infeasible cases
            end
        end
    end
end
function V_vax = utility_vax(I_grid, R_grid, PAR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UTILITY_VAX  Compute the continuation value of welfare after arrival
    % of the vaccine. Called by simulations.m, backward_induct_comm.m, and
    % backward_induct_nocomm.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   I_grid = initial share of infected individuals.
    %   R_grid = initial share of recovered individuals.
    %   PAR = structure array of parameters.
    %
    % Outputs:
    %   V_vax = continuation value after arrival of the vaccine.
    
    
    %%% initialization
    % setup grids of health state and lockdown policy
    [Ig, Rg] = meshgrid(I_grid, R_grid); % set up grid with I and R. Here the dim of grid is n_R*n_I

    [~, I_prime, R_prime, ~] = sir_model(Ig, Rg, 1, PAR, true); % next period number of each variety of people under optimal policy of no lockdown once the vaccine has arrived

    % functional forms
    A = PAR.TFP1;
    S = @(I, R) 1 - I - R*(1 + PAR.RHO4/PAR.RHO3); % share susceptible (S) from rearranging SIR model
    l_star = @(I, R) S(I, R) + PAR.GAMMA*I + R; % effective labor (\ell) incl. lockdown (1 - L)
    x_star = @(I, R) (PAR.ALPHA.*A/PAR.INTEREST_W).^(1/(1 - PAR.ALPHA)).*l_star(I, R); % investment (x) from no-arbitrage condition
    w_star = @(I, R) (1 - PAR.ALPHA)*A*(x_star(I, R)./l_star(I, R)).^PAR.ALPHA;
    c_star = @(I, R) w_star(I, R).*l_star(I, R)./(S(I, R) + I + R);
    u_star = @(I, R) (S(I, R) + I + R).*(log(c_star(I, R)) + PAR.NU);

    
    %%% value function iteration
    % set initial value, assume to be SS value
    V_vax = u_star(Ig, Rg)/(1 - PAR.BETA); % set initial value, assume to be SS value

    n_I = length(I_grid);
    n_R = length(R_grid);

    vinterp = u_star(Ig(:), Rg(:)) + PAR.BETA*interp2(I_grid, R_grid, V_vax, I_prime(:), R_prime(:), 'linear'); 
    vinterp = reshape(vinterp, n_R, n_I); % update

    crit = 1e-8; 

    % iteration
    while max(max(abs(vinterp - V_vax))) >= crit
        V_vax = vinterp;
        vinterp = u_star(Ig(:), Rg(:)) + PAR.BETA*interp2(I_grid, R_grid, V_vax, I_prime(:), R_prime(:), 'linear'); % update
        vinterp = reshape(vinterp, n_R, n_I);
    end
end
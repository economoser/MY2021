function [Vold, pol, Vstore, Pstore] = backward_induct_comm(I_grid, R_grid, L_grid, tt, TT, welfare_function, PAR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARD_INDUCT_COMM  Return the continuation value and the optimal
    % lockdown policy from time tt to TT for a government with commitment.
    % Called by simulations.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   I_grid = share of infected individuals.
    %   R_grid = share of recovered individuals.
    %   L_grid = lockdown policy.
    %   tt = initial period.
    %   TT = terminal period.
    %   welfare_function = welfare in each period with commitment.
    %   PAR = structure array of parameters.
    %
    % Outputs:
    %   Vold = continuation value at beginning given health state.
    %   pol = optimal lockdown policy at beginning given health state.
    %   Vstore = continuation value in each period given health state.
    %   Pstore = optimal lockdown policy in each period given health state.

    
    %%% initialization
    % setup grids of health state and lockdown policy
    n_I = length(I_grid);
    n_R = length(R_grid);
    n_L = length(L_grid);
    [Ig, Rg, Lg] = meshgrid(I_grid, R_grid, L_grid); % Note: The three LHS objects are each of dimension n_R x n_I x n_L. The MATLAB command -interp2- uses the same convention.

    % set initial value for value and optimal policy
    Vstore = NaN(n_I, n_R, TT - tt + 1);
    Pstore = NaN(n_I, n_R, TT - tt + 1);
    
    % calculate continuation value at time TT, when vaccine arrives.
    Vold = utility_vax(I_grid, R_grid, PAR);

    Vstore(:, :, end) = Vold';
    Pstore(:, :, end) = ones(n_I, n_R); % after vax no quarantine.
    
    
    %%% backward induction
    for time = TT:-1:tt
        [~, I_prime, R_prime, ~] = sir_model(Ig, Rg, Lg, PAR, false);
        
        vinterp = welfare_function(Ig(:), Rg(:), Lg(:)) + PAR.BETA*interp2(I_grid, R_grid, Vold, I_prime(:), R_prime(:), 'linear'); % current utility + continuation
        vinterp = reshape(vinterp, n_R, n_I, n_L);
        
        [Vold, pol] = max(vinterp, [], 3); % optimal lockdown policy and corresponding continuation of welfare

        Vstore(:, :, time - tt + 1) = Vold';
        Pstore(:, :, time - tt + 1) = L_grid(pol');
    end
end
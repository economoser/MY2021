function [Vold, pol, Vstore, Pstore] = backward_induct_nocomm(I_grid, R_grid, L_grid, tt, TT, welfare_function_no_commit, welfare_function_commit, PAR)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BACKWARD_INDUCT_NOCOMM  Return the continuation value and the optimal
    % lockdown policy from time tt to TT for a government without
    % commitment. Called by simulations.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   I_grid = share of infected individuals.
    %   R_grid = share of recovered individuals.
    %   L_grid = lockdown policy.
    %   tt = initial period.
    %   TT = terminal period.
    %   welfare_function_no_commit = welfare in each period with the
    %      government lack of commitment, inputs are the health state and
    %      lockdown policy.
    %   welfare_function_commit = welfare in each period with commitment, 
    %      inputs are the health state and lockdown policy.
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

    [Ig2, Rg2] = meshgrid(I_grid, R_grid);
    
    % set initial value for value and optimal policy
    Vstore = NaN(n_I, n_R, TT - tt + 1);
    Pstore = NaN(n_I, n_R, TT - tt + 1);
    
    % calculate continuation value at time TT, when vaccine arrives.
    Vold = utility_vax(I_grid, R_grid, PAR); % continuation of welfare after vax

    Vstore(:, :, end) = Vold';
    Pstore(:, :, end) = ones(n_I,n_R); % after vax no quarantine.
    
    
    %%% backward induction
    for time = TT:-1:tt
        [~, I_prime, R_prime, ~] = sir_model(Ig, Rg, Lg, PAR, false);
        
        vinterp = welfare_function_no_commit(Ig(:), Rg(:), Lg(:)) + PAR.BETA*interp2(I_grid, R_grid, Vold, I_prime(:), R_prime(:), 'linear'); % current utility + continuation
        vinterp = reshape(vinterp, n_R, n_I, n_L); 
        
        [Vold, pol] = max(vinterp, [], 3); % optimal lockdown policy under lack of commitment and corresponding continuation of welfare

        [~, I_prime, R_prime, ~] = sir_model(Ig2, Rg2, L_grid(pol), PAR, false);
        
        pol = L_grid(pol);
        vinterp = welfare_function_commit(Ig2(:), Rg2(:), pol(:)) + PAR.BETA*interp2(I_grid, R_grid, Vold, I_prime(:), R_prime(:), 'linear');
        Vold = reshape(vinterp, n_R, n_I); 

        Vstore(:, :, time - tt + 1) = Vold';
        Pstore(:, :, time - tt + 1) = pol';
    end  
end
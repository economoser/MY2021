function [lockdown_disc, lockdown_disc_pand, lockdown_disc_year, lockdown_thresh, lockdown_thresh_pand, lockdown_thresh_year, gdp, gdp_pand, gdp_year, deaths, deaths_pand, deaths_year, alive_disc, alive_disc_pand, alive_disc_year, V, V_pand, V_year, cons_equiv_loss, cons_equiv_loss_pand, cons_equiv_loss_year] = simulations(PAR, DIR_SIM, replace_name, replace_value) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATIONS  Simulate the evolution of key variables and return
    % summary statistics. Called by master.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   PAR = structure array of parameters.
    %   DIR_SIM = directory in which to store simulation results.
    %   replace_name = parameters to be adjusted for comparative statics.
    %   replace_value = new values of the parameters to be adjusted.
    %
    % Outputs:
    %   lockdown_disc = NPV of lockdown policy in infinite-horizon economy
    %      (t = 1, 2, ..., +infinity).
    %   lockdown_disc_pand = NPV of lockdown policy during pandemic period
    %      (t = 1, 2, ..., TTP).
    %   lockdown_disc_year = NPV of lockdown policy during first year (t =
    %      1, 2, ..., 52).
    %   lockdown_thresh = NPV of lockdown policy in infinite-horizon
    %      economy (t = 1, 2, ..., +infinity).
    %   lockdown_thresh_pand = NPV of lockdown policy during pandemic
    %      period (t = 1, 2, ..., TTP).
    %   lockdown_thresh_year = NPV of lockdown policy during first year (t
    %      = 1, 2, ..., 52).
    %   gdp = NPV of GDP in infinite-horizon economy (t = 1, 2, ...,
    %      +infinity).
    %   gdp_pand = NPV of GDP during pandemic period (t = 1, 2, ..., TTP).
    %   gdp_year = NPV of GDP during first year (t = 1, 2, ..., 52).
    %   deaths = share of population who die in infinite-horizon economy (t
    %      = 1, 2, ..., +infinity).
    %   deaths_pand = share of population who die during pandemic period (t
    %      = 1, 2, ..., TTP).
    %   deaths_year = share of population who die during first year (t = 1,
    %      2, ..., 52).
    %   alive_disc = share of population who die in infinite-horizon
    %      economy (t = 1, 2, ..., +infinity).
    %   alive_disc_pand = share of population who die during pandemic
    %      period (t = 1, 2, ..., TTP).
    %   alive_disc_year = share of population who die during first year (t
    %      = 1, 2, ..., 52).
    %   V = value in infinite-horizon economy (t = 1, 2, ..., +infinity).
    %   V_pand = value during pandemic period (t = 1, 2, ..., TTP).
    %   V_year = value during first year (t = 1, 2, ..., 52).
    %   cons_equiv_loss = consumption-equivalent welfare losses from the
    %      pandemic in infinite-horizon economy (t = 1, 2, ..., +infinity).
    %   cons_equiv_loss_pand = consumption-equivalent welfare losses from
    %      the pandemic during pandemic period (t = 1, 2, ..., TTP).
    %   cons_equiv_loss_year = consumption-equivalent welfare losses from
    %      the pandemic during first year (t = 1, 2, ..., 52).

    
    %%% set parameter values
    TT = PAR.TT;
    TTP = PAR.TTP;
    
    
    %%% replace parameter values of requested for comparative statics
    if ~isempty(replace_name)
        PAR.(replace_name) = replace_value;
        suffix_out = sprintf(['_', replace_name, '_%10.8f'], replace_value);
        suffix_out = strrep(suffix_out, '.', 'p');
    else
        suffix_out = '_baseline';
    end
    
    
    %%% functional forms
    % technology
    A = @(I) PAR.TFP1*(1 - PAR.TFP2*I.^PAR.TFP3); % TFP (A) from functional form assumption
    S = @(I, R) 1 - I - R*(1 + PAR.RHO4/PAR.RHO3); % share susceptible (S) from rearranging SIR model
    l_star = @(I, R, L) L.*(S(I, R) + PAR.GAMMA*I + R); % effective labor (\ell) incl. lockdown (1 - L)
    x_star = @(I, R, L) (PAR.ALPHA.*A(I)/PAR.INTEREST_W).^(1/(1 - PAR.ALPHA)).*l_star(I, R, L); % investment (x) from no-arbitrage condition
    w_star = @(I, R, L) (1 - PAR.ALPHA)*A(I).*(x_star(I, R, L)./l_star(I, R, L)).^PAR.ALPHA; % equilibrium wage (w) from competitive-labor-market condition
    c_star = @(I, R, L) w_star(I, R, L).*l_star(I, R, L)./(S(I, R) + I + R); % per capita consumption (c/(S + I + R)) from budget constraint
    y_star = @(I, R, L) A(I).*x_star(I, R, L).^(PAR.ALPHA).*l_star(I, R, L).^(1 - PAR.ALPHA); % output (y) from functional form assumption
    
    % period utility
    if PAR.SIGMA == 1.0
      u_star = @(I, R, L) log(c_star(I, R, L)) + PAR.NU; % total flow utility
      u_g = @(I, R, L) (1 - PAR.ALPHA)*log(c_star(I, R, L)) + PAR.NU; % utility with no commitment
      u_e = @(I, R, L) log(c_star(I, R, L)) + PAR.NU; % utility with commitment
    else
      u_star = @(I, R, L) (c_star(I, R, L).^(1 - PAR.SIGMA) - 1)./(1 - PAR.SIGMA) + PAR.NU; % period utility
      u_g = @(I, R, L) (1 - PAR.ALPHA)*(c_star(I, R, L).^(1 - PAR.SIGMA) - 1)./(1 - PAR.SIGMA) + PAR.NU; % utility with no commitment
      u_e = @(I, R, L) (c_star(I, R, L).^(1 - PAR.SIGMA) - 1)./(1 - PAR.SIGMA) + PAR.NU; % utility with commitment
    end
    
    % period welfare 
    welfare_star = @(I, R, L) (S(I, R) + I + R).*u_star(I, R, L);
    welfare_g = @(I, R, L) (S(I, R) + I + R).*u_g(I, R, L);
    welfare_e = @(I, R, L) (S(I, R) + I + R).*u_e(I, R, L);
    
    % steady-state per-capita consumption without a pandemic
    c_tilde = (1 - PAR.ALPHA)*A(0)^(1/(1 - PAR.ALPHA))*(PAR.ALPHA/PAR.INTEREST_W)^(PAR.ALPHA/(1 - PAR.ALPHA));
    
    
    %%% initialize vectors to store simulation results
    % no pandemic
    SH = NaN(TTP, 1);
    IH = NaN(TTP, 1);
    RH = NaN(TTP, 1);
    DH = NaN(TTP, 1);
    LH = NaN(TTP, 1);
    
    % commitment
    SC = NaN(TTP, 1);
    IC = NaN(TTP, 1);
    RC = NaN(TTP, 1);
    DC = NaN(TTP, 1);
    LC = NaN(TTP, 1);
    
    % no commitment
    SN = NaN(TTP, 1);
    IN = NaN(TTP, 1);
    LN = NaN(TTP, 1);
    RN = NaN(TTP, 1);
    DN = NaN(TTP, 1);
    
    % no lockdown
    S0 = NaN(TTP, 1);
    I0 = NaN(TTP, 1);
    R0 = NaN(TTP, 1);
    D0 = NaN(TTP, 1);
    L0 = NaN(TTP, 1);

    
    %%% grids
    R_grid_max = PAR.RHO3/(PAR.RHO3 + PAR.RHO4); % R has a maximal value as a function of rho3 and rho4
    I_grid = linspace(0, 1, PAR.N_GRID_I);
    R_grid = linspace(0, R_grid_max, PAR.N_GRID_R);
    L_grid = linspace(0, 1, PAR.N_GRID_L); % This is 1 - lockdown policy [0: nobody works, 1: no lockdown]
    I_grid(1) = 0;
    R_grid(1) = 0;
    L_grid(1) = L_grid(2)/2;
    I_grid(end) = 1;
    R_grid(end) = R_grid_max;
    L_grid(end) = 1;

    I_grid = I_grid.^PAR.GRID_POWER;
    R_grid = (R_grid/R_grid_max).^PAR.GRID_POWER*R_grid_max;
    L_grid = L_grid.^PAR.GRID_POWER;
    
    
    %%% find problem solutions
    % backward induction
    [~, ~, ~, PC] = backward_induct_comm(I_grid, R_grid, L_grid, 1, TT, welfare_star, PAR); % solve problem with commitment
    [~, ~, ~, PN] = backward_induct_nocomm(I_grid, R_grid, L_grid, 1, TT, welfare_g, welfare_e, PAR); % solve problem with no commitment
    
    % store parameter values
    RHO1_store = PAR.RHO1;
    RHO2_store = PAR.RHO2;
    
    
    %%% simulate economy for t = 1, ..., TTP
    for tt = 1:TTP
      if tt == 1
        ICtt = PAR.I_INIT;
        RCtt = 0;
        INtt = PAR.I_INIT;
        RNtt = 0;
        I0tt = PAR.I_INIT;
        R0tt = 0;
      else
        ICtt = IC(tt - 1);
        RCtt = RC(tt - 1);
        INtt = IN(tt - 1);
        RNtt = RN(tt - 1);
        I0tt = I0(tt - 1);
        R0tt = R0(tt - 1);
      end  
      
      % no pandemic
      SH(tt) = 1;
      IH(tt) = 0;
      RH(tt) = 0;
      DH(tt) = 0;
      LH(tt) = 1;
      
      % commitment
      if tt > TT
        PAR.RHO1 = 0;
        PAR.RHO2 = 0;
        LC(tt) = 1;
        [SC(tt), IC(tt), RC(tt), DC(tt)] = sir_model(ICtt, RCtt, LC(tt), PAR, true);
      else
        LC(tt) = interp2(I_grid, R_grid, PC(:, :, tt)', ICtt, RCtt, 'linear');
        LC(tt) = min(max(LC(tt), 0), 1); % constrain to be in [0,1]
        [SC(tt), IC(tt), RC(tt), DC(tt)] = sir_model(ICtt, RCtt, LC(tt), PAR, false);
      end
      
      % no commitment
      if tt > TT
        LN(tt) = 1;
        [SN(tt), IN(tt), RN(tt), DN(tt)] = sir_model(INtt, RNtt, LN(tt), PAR, true);
      else
        LN(tt) = interp2(I_grid, R_grid, PN(:, :, tt)', INtt, RNtt, 'linear');
        LN(tt) = min(max(LN(tt), 0), 1); % constrain to be in [0,1]
        [SN(tt), IN(tt), RN(tt), DN(tt)] = sir_model(INtt, RNtt, LN(tt), PAR, false);
      end
      
      % no lockdown
      L0(tt) = 1.0;
      if tt>TT
        [S0(tt), I0(tt), R0(tt), D0(tt)] = sir_model(I0tt, R0tt, L0(tt), PAR, true);
      else
        [S0(tt), I0(tt), R0(tt), D0(tt)] = sir_model(I0tt, R0tt, L0(tt), PAR, false);
      end
    end
    
    % steady-state health state with no pandemic
    DH_ss = 0;
    
    % steady-state health state with commitment
    LC_ss = 1;
    SC_ss = SC(end);
    IC_ss = IC(end);
    RC_ss = RC(end);
    DC_ss = DC(end);
    diff = max([SC_ss - SC(end - 1), IC_ss - IC(end - 1), RC_ss - RC(end - 1), DC_ss - DC(end - 1)]);
    while diff > 10e-6
        SC_ss_old = SC_ss;
        IC_ss_old = IC_ss;
        RC_ss_old = RC_ss;
        DC_ss_old = DC_ss;
        [SC_ss, IC_ss, RC_ss, DC_ss] = sir_model(IC_ss, RC_ss, LC_ss, PAR, true);
        diff = max([SC_ss - SC_ss_old, IC_ss - IC_ss_old, RC_ss - RC_ss_old, DC_ss - DC_ss_old]);
    end
    
    % steady-state health state with no commitment
    LN_ss = 1;
    SN_ss = SN(end);
    IN_ss = IN(end);
    RN_ss = RN(end);
    DN_ss = DN(end);
    diff = max([SN_ss - SN(end - 1), IN_ss - IN(end - 1), RN_ss - RN(end - 1), DN_ss - DN(end - 1)]);
    while diff > 10e-6
        SN_ss_old = SN_ss;
        IN_ss_old = IN_ss;
        RN_ss_old = RN_ss;
        DN_ss_old = DN_ss;
        [SN_ss, IN_ss, RN_ss, DN_ss] = sir_model(IN_ss, RN_ss, LN_ss, PAR, true);
        diff = max([SN_ss - SN_ss_old, IN_ss - IN_ss_old, RN_ss - RN_ss_old, DN_ss - DN_ss_old]);
    end
    
    % steady-state health state with no lockdown
    L0_ss = 1;
    S0_ss = S0(end);
    I0_ss = I0(end);
    R0_ss = R0(end);
    D0_ss = D0(end);
    diff = max([S0_ss - S0(end - 1), I0_ss - I0(end - 1), R0_ss - R0(end - 1), D0_ss - D0(end - 1)]);
    while diff > 10e-6
        S0_ss_old = S0_ss;
        I0_ss_old = I0_ss;
        R0_ss_old = R0_ss;
        D0_ss_old = D0_ss;
        [S0_ss, I0_ss, R0_ss, D0_ss] = sir_model(I0_ss, R0_ss, L0_ss, PAR, true);
        diff = max([S0_ss - S0_ss_old, I0_ss - I0_ss_old, R0_ss - R0_ss_old, D0_ss - D0_ss_old]);
    end
    
    % recover stored parameter values
    PAR.RHO1 = RHO1_store;
    PAR.RHO2 = RHO2_store;

    
    %%% table with results from simulation of time series
    if strcmp(suffix_out,'_baseline')
        results_time_series = table( ...
            (1:TTP)', (1:TTP)' > TT, ...
            1 - L0, 1 - LC, 1 - LN, ...
            y_star(I0, R0, L0), y_star(IC, RC, LC), y_star(IN, RN, LN), ...
            x_star(I0, R0, L0), x_star(IC, RC, LC), x_star(IN, RN, LN), ...
            log(c_star(I0, R0, L0)), log(c_star(IC, RC, LC)), log(c_star(IN, RN, LN)), ...
            repmat(PAR.NU, size(c_star(I0, R0, L0))), repmat(PAR.NU, size(c_star(I0, R0, L0))), repmat(PAR.NU, size(c_star(I0, R0, L0))), ...
            S0, SC, SN, ...
            I0, IC, IN, ...
            R0, RC, RN, ...
            D0, DC, DN ...
        );
        results_time_series.Properties.VariableNames = { ...
            'Time', 'Vaccine', ...
            'lock0','lockC','lockN', ...
            'Y0','YC','YN', ...
            'X0','XC','XN', ...
            'LOGC_0','LOGC_C','LOGC_N', ...
            'NU_0','NU_C','NU_N', ...
            'S0','SC','SN', ...
            'I0','IC','IN', ...
            'R0','RC','RN', ...
            'D0','DC','DN' ...
        };
        FILE_TIME_SERIES = [DIR_SIM, '/sim_time_series.xlsx'];
        if exist(FILE_TIME_SERIES, 'file')
            delete(FILE_TIME_SERIES)
        end
        writetable(results_time_series, FILE_TIME_SERIES);
    end
    
    
    %%% summary statistics
    disc_factor = 1/(1 + PAR.INTEREST_W);
    y_vax = gdp_vax(I_grid, R_grid, PAR);
    
    lockdown_disc = struct(); % NPV of lockdown policy in infinite-horizon economy (t = 1, 2, ..., +infinity)
    lockdown_disc.('no_pand') = sum(disc_factor.^(0:TTP - 1).*(1 - LH'));
    lockdown_disc.('commit') = sum(disc_factor.^(0:TTP - 1).*(1 - LC'));
    lockdown_disc.('no_commit') = sum(disc_factor.^(0:TTP - 1).*(1 - LN'));
    lockdown_disc.('no_lock') = sum(disc_factor.^(0:TTP - 1).*(1 - L0'));

    lockdown_disc_pand = struct(); % NPV of lockdown policy during pandemic period (t = 1, 2, ..., TTP)
    lockdown_disc_pand.('no_pand') = sum(disc_factor.^(0:TTP - 1).*(1 - LH'));
    lockdown_disc_pand.('commit') = sum(disc_factor.^(0:TTP - 1).*(1 - LC'));
    lockdown_disc_pand.('no_commit') = sum(disc_factor.^(0:TTP - 1).*(1 - LN'));
    lockdown_disc_pand.('no_lock') = sum(disc_factor.^(0:TTP - 1).*(1 - L0'));
    
    lockdown_disc_year = struct(); % NPV of lockdown policy during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        lockdown_disc_year.('no_pand') = sum(disc_factor.^(0:51).*(1 - LH(1:52)'));
        lockdown_disc_year.('commit') = sum(disc_factor.^(0:51).*(1 - LC(1:52)'));
        lockdown_disc_year.('no_commit') = sum(disc_factor.^(0:51).*(1 - LN(1:52)'));
        lockdown_disc_year.('no_lock') = sum(disc_factor.^(0:51).*(1 - L0(1:52)'));
    else
        lockdown_disc_year.('no_pand') = 0;
        lockdown_disc_year.('commit') = 0;
        lockdown_disc_year.('no_commit') = 0;
        lockdown_disc_year.('no_lock') = 0; 
    end
    
    lockdown_thresh = struct(); % NPV of lockdown policy in infinite-horizon economy (t = 1, 2, ..., +infinity)
    thresh = .05; % minimum lockdown threshold above which to count number of days in lockdown
    lockdown_thresh.('no_pand') = sum(1 - LH' > thresh);
    lockdown_thresh.('commit') = sum(1 - LC' > thresh);
    lockdown_thresh.('no_commit') = sum(1 - LN' > thresh);
    lockdown_thresh.('no_lock') = sum(1 - L0' > thresh);

    lockdown_thresh_pand = struct(); % NPV of lockdown policy during pandemic period (t = 1, 2, ..., TTP)
    lockdown_thresh_pand.('no_pand') = sum(1 - LH' > thresh);
    lockdown_thresh_pand.('commit') = sum(1 - LC' > thresh);
    lockdown_thresh_pand.('no_commit') = sum(1 - LN' > thresh);
    lockdown_thresh_pand.('no_lock') = sum(1 - L0' > thresh);
    
    lockdown_thresh_year = struct(); % NPV of lockdown policy during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        lockdown_thresh_year.('no_pand') = sum(1 - LH(1:52)' > thresh);
        lockdown_thresh_year.('commit') = sum(1 - LC(1:52)' > thresh);
        lockdown_thresh_year.('no_commit') = sum(1 - LN(1:52)' > thresh);
        lockdown_thresh_year.('no_lock') = sum(1 - L0(1:52)' > thresh);
    else
        lockdown_thresh_year.('no_pand') = 0;
        lockdown_thresh_year.('commit') = 0;
        lockdown_thresh_year.('no_commit') = 0;
        lockdown_thresh_year.('no_lock') = 0; 
    end
    
    gdp = struct(); % NPV of GDP in infinite-horizon economy (t = 1, 2, ..., +infinity)
    gdp.('no_pand') = sum(disc_factor.^(0:TTP - 1).*y_star(IH, RH, LH)') + disc_factor^TTP*interp2(I_grid, R_grid, y_vax, IH(end), RH(end), 'linear');
    gdp.('commit') = sum(disc_factor.^(0:TTP - 1).*y_star(IC, RC, LC)') + disc_factor^TTP*interp2(I_grid, R_grid, y_vax, IC(end), RC(end), 'linear');
    gdp.('no_commit') = sum(disc_factor.^(0:TTP - 1).*y_star(IN, RN, LN)') + disc_factor^TTP*interp2(I_grid, R_grid, y_vax, IN(end), RN(end), 'linear');
    gdp.('no_lock') = sum(disc_factor.^(0:TTP - 1).*y_star(I0, R0, L0)') + disc_factor^TTP*interp2(I_grid, R_grid, y_vax, I0(end), R0(end), 'linear');

    gdp_pand = struct(); % NPV of GDP during pandemic period (t = 1, 2, ..., TTP)
    gdp_pand.('no_pand') = sum(disc_factor.^(0:TTP - 1).*y_star(IH, RH, LH)');
    gdp_pand.('commit') = sum(disc_factor.^(0:TTP - 1).*y_star(IC, RC, LC)');
    gdp_pand.('no_commit') = sum(disc_factor.^(0:TTP - 1).*y_star(IN, RN, LN)');
    gdp_pand.('no_lock') = sum(disc_factor.^(0:TTP - 1).*y_star(I0, R0, L0)');
    
    gdp_year = struct(); % NPV of GDP during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        gdp_year.('no_pand') = sum(disc_factor.^(0:51).*y_star(IH(1:52), RH(1:52), LH(1:52))');
        gdp_year.('commit') = sum(disc_factor.^(0:51).*y_star(IC(1:52), RC(1:52), LC(1:52))');
        gdp_year.('no_commit') = sum(disc_factor.^(0:51).*y_star(IN(1:52), RN(1:52), LN(1:52))');
        gdp_year.('no_lock') = sum(disc_factor.^(0:51).*y_star(I0(1:52), R0(1:52), L0(1:52))');
    else
        gdp_year.('no_pand') = 0;
        gdp_year.('commit') = 0;
        gdp_year.('no_commit') = 0;
        gdp_year.('no_lock') = 0; 
    end

    deaths = struct(); % share of population who die in infinite-horizon economy (t = 1, 2, ..., +infinity)
    deaths.('no_pand') = DH_ss;
    deaths.('commit') = DC_ss;
    deaths.('no_commit') = DN_ss;
    deaths.('no_lock') = D0_ss;
    
    deaths_pand = struct(); % share of population who die during pandemic period (t = 1, 2, ..., TTP)
    deaths_pand.('no_pand') = DH(TTP);
    deaths_pand.('commit') = DC(TTP);
    deaths_pand.('no_commit') = DN(TTP);
    deaths_pand.('no_lock') = D0(TTP);
    
    deaths_year = struct(); % share of population who die during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        deaths_year.('no_pand') = DH(52);
        deaths_year.('commit') = DC(52);
        deaths_year.('no_commit') = DN(52);
        deaths_year.('no_lock') = D0(52);
    else
        deaths_year.('no_pand') = 0;
        deaths_year.('commit') = 0;
        deaths_year.('no_commit') = 0;
        deaths_year.('no_lock') = 0;
    end
    
    alive_disc = struct(); % share of population who die in infinite-horizon economy (t = 1, 2, ..., +infinity)
    alive_disc.('no_pand') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DH')) + PAR.BETA^TTP*(1 - DH_ss)/(1 - PAR.BETA);
    alive_disc.('commit') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DC')) + PAR.BETA^TTP*(1 - DC_ss)/(1 - PAR.BETA);
    alive_disc.('no_commit') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DN')) + PAR.BETA^TTP*(1 - DN_ss)/(1 - PAR.BETA);
    alive_disc.('no_lock') = sum(PAR.BETA.^(0:TTP - 1).*(1 - D0')) + PAR.BETA^TTP*(1 - D0_ss)/(1 - PAR.BETA);
    
    alive_disc_pand = struct(); % share of population who die during pandemic period (t = 1, 2, ..., TTP)
    alive_disc_pand.('no_pand') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DH(1:TTP)'));
    alive_disc_pand.('commit') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DC(1:TTP)'));
    alive_disc_pand.('no_commit') = sum(PAR.BETA.^(0:TTP - 1).*(1 - DN(1:TTP)'));
    alive_disc_pand.('no_lock') = sum(PAR.BETA.^(0:TTP - 1).*(1 - D0(1:TTP)'));
    
    alive_disc_year = struct(); % share of population who die during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        alive_disc_year.('no_pand') = sum(PAR.BETA.^(0:51).*(1 - DH(1:52)'));
        alive_disc_year.('commit') = sum(PAR.BETA.^(0:51).*(1 - DC(1:52)'));
        alive_disc_year.('no_commit') = sum(PAR.BETA.^(0:51).*(1 - DN(1:52)'));
        alive_disc_year.('no_lock') = sum(PAR.BETA.^(0:51).*(1 - D0(1:52)'));
    else
        alive_disc_year.('no_pand') = 0;
        alive_disc_year.('commit') = 0;
        alive_disc_year.('no_commit') = 0;
        alive_disc_year.('no_lock') = 0;
    end
    
    Vold = utility_vax(I_grid, R_grid, PAR);
    
    V = struct(); % value in infinite-horizon economy (t = 1, 2, ..., +infinity)
    V.('no_pand') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IH, RH, LH)') + PAR.BETA^TTP*interp2(I_grid, R_grid, Vold, IH(end),RH(end), 'linear');
    V.('commit') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IC, RC, LC)') + PAR.BETA^TTP*interp2(I_grid, R_grid, Vold, IC(end),RC(end), 'linear');
    V.('no_commit') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IN, RN, LN)') + PAR.BETA^TTP*interp2(I_grid, R_grid, Vold, IN(end),RN(end), 'linear');
    V.('no_lock') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(I0, R0, L0)') + PAR.BETA^TTP*interp2(I_grid, R_grid, Vold, I0(end),R0(end), 'linear');

    V_pand = struct(); % value during pandemic period (t = 1, 2, ..., TTP)
    V_pand.('no_pand') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IH, RH, LH)');
    V_pand.('commit') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IC, RC, LC)');
    V_pand.('no_commit') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(IN, RN, LN)');
    V_pand.('no_lock') = sum(PAR.BETA.^(0:TTP - 1).*welfare_star(I0, R0, L0)');
    
    V_year = struct(); % value during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        V_year.('no_pand') = sum(PAR.BETA.^(0:51).*welfare_star(IH(1:52), RH(1:52), LH(1:52))');
        V_year.('commit') = sum(PAR.BETA.^(0:51).*welfare_star(IC(1:52), RC(1:52), LC(1:52))');
        V_year.('no_commit') = sum(PAR.BETA.^(0:51).*welfare_star(IN(1:52), RN(1:52), LN(1:52))');
        V_year.('no_lock') = sum(PAR.BETA.^(0:51).*welfare_star(I0(1:52), R0(1:52), L0(1:52))');
    else
        V_year.('no_pand') = 0;
        V_year.('commit') = 0;
        V_year.('no_commit') = 0;
        V_year.('no_lock') = 0; 
    end
    
    cons_equiv_loss = struct(); % consumption-equivalent welfare losses from the pandemic in infinite-horizon economy (t = 1, 2, ..., +infinity)
    if PAR.SIGMA == 1.0
        cons_equiv_loss.('no_pand') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_pand') - PAR.NU)/c_tilde;
        cons_equiv_loss.('commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('commit') - PAR.NU)/c_tilde;
        cons_equiv_loss.('no_commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_commit') - PAR.NU)/c_tilde;
        cons_equiv_loss.('no_lock') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_lock') - PAR.NU)/c_tilde;
    else
        cons_equiv_loss.('no_pand')   = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_pand') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss.('commit')    = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss.('no_commit') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss.('no_lock')   = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^inf)*V.('no_lock') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
    end
    
    cons_equiv_loss_pand = struct(); % consumption-equivalent welfare losses from the pandemic during pandemic period (t = 1, 2, ..., TTP)
    if PAR.SIGMA == 1.0
        cons_equiv_loss_pand.('no_pand') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_pand') - PAR.NU)/c_tilde;
        cons_equiv_loss_pand.('commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('commit') - PAR.NU)/c_tilde;
        cons_equiv_loss_pand.('no_commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_commit') - PAR.NU)/c_tilde;
        cons_equiv_loss_pand.('no_lock') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_lock') - PAR.NU)/c_tilde;
    else
        cons_equiv_loss_pand.('no_pand') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_pand') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss_pand.('commit') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss_pand.('no_commit') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        cons_equiv_loss_pand.('no_lock') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^TTP)*V_pand.('no_lock') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
    end
    
    cons_equiv_loss_year = struct(); % consumption-equivalent welfare losses from the pandemic during first year (t = 1, 2, ..., 52)
    if TTP >= 52
        if PAR.SIGMA == 1.0
            cons_equiv_loss_year.('no_pand') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_pand') - PAR.NU)/c_tilde;
            cons_equiv_loss_year.('commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('commit') - PAR.NU)/c_tilde;
            cons_equiv_loss_year.('no_commit') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_commit') - PAR.NU)/c_tilde;
            cons_equiv_loss_year.('no_lock') = 1 - exp((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_lock') - PAR.NU)/c_tilde;
        else
            cons_equiv_loss_year.('no_pand') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_pand') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
            cons_equiv_loss_year.('commit') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
            cons_equiv_loss_year.('no_commit') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_commit') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
            cons_equiv_loss_year.('no_lock') = 1 - ((1 - PAR.SIGMA)*((1 - PAR.BETA)/(1 - PAR.BETA^52)*V_year.('no_lock') - PAR.NU) + 1)^(1/(1 - PAR.SIGMA))/c_tilde;
        end
    else
        cons_equiv_loss_year.('no_pand') = 0;
        cons_equiv_loss_year.('commit') = 0;
        cons_equiv_loss_year.('no_commit') = 0;
        cons_equiv_loss_year.('no_lock') = 0;
    end
end
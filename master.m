%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASTER  Simulate an economy subject to optimal lockdown policy with and
% without commitment for Moser and Yared (2021). MATLAB master file that
% can also be called from Stata in master.do -- run this first!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% initial housekeeping
% clear memory and display
clear
clc

% start timer
tic

% section switches
switch_sim = 1; % 0 = do not run baseline simulation, 1 = run baseline simulation
switch_comp_binary = 1; % 0 = do not run low-high comparative statics, 1 = run low-high comparative statics
switch_comp_multiple = 1; % 0 = do not run multiple comparative statics, 1 = run multiple comparative statics

% set main directory
dir_list = { ...
    'SET_MAIN_DIRECTORY_HERE' ...
};
for n = 1:length(dir_list)
    dir = dir_list{n};
    if exist(dir, 'dir')
       DIR_MAIN = dir;
    end
end
if ~exist('DIR_MAIN', 'var')
   error('USER ERROR: Must specify main directory in line 23!')
end

% create subdirectories
DIR_CODE = [DIR_MAIN, '/_code'];
if ~exist(DIR_CODE, 'dir')
    mkdir(DIR_CODE)
end
DIR_LOG = [DIR_MAIN, '/_log'];
if ~exist(DIR_LOG, 'dir')
    mkdir(DIR_LOG)
end
DIR_RESULTS = [DIR_MAIN, '/_results'];
if ~exist(DIR_RESULTS, 'dir')
    mkdir(DIR_RESULTS)
end
DIR_FIG = [DIR_RESULTS, '/_figures'];
if ~exist(DIR_FIG, 'dir')
    mkdir(DIR_FIG)
end
DIR_PAR = [DIR_RESULTS, '/_parallelization'];
if exist(DIR_PAR, 'dir')
    rmdir(DIR_PAR, 's')
end
mkdir(DIR_PAR)
DIR_SIM = [DIR_RESULTS, '/_simulations'];
if ~exist(DIR_SIM, 'dir')
    mkdir(DIR_SIM)
end

% set working directory
cd(DIR_CODE)

% open log file
FILE_LOG = [DIR_LOG, '/log_covid_matlab.log'];
if exist(FILE_LOG, 'file')
    delete(FILE_LOG)
end
diary(FILE_LOG)

% start parallel pool
delete(gcp('nocreate'))
c = parcluster('local');
c.JobStorageLocation = DIR_PAR; % recommended to avoid storing conflicting job information (see section "Running Multiple PCT MATLAB Jobs" at https://rcc.uchicago.edu/docs/software/environments/matlab/)
nw = c.NumWorkers; % number of available cores
parpool(c, nw, 'IdleTimeout', Inf); % assign all cores to MATLAB

% baseline and binary comparative statics
l_min = 2 - max(switch_sim, switch_comp_multiple);
l_max = switch_comp_binary*16 + (1 - switch_comp_binary)*1;
for l = l_min:l_max
    if switch_comp_binary && l < 16
        fprintf(['\n\n\n\n\n\n\n\n\n\n*** loop number = ', num2str(l), '\n'])
    end
    
    
    %%% model parameters
    % calibration targets
    length_infection = 14/7; % mean length of infection (weeks = days/(days/week)) from https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html -- default = 2.0
    mortality_rate = 0.058204761; % mortality rate conditional on infection (deaths/case) from https://covidtracking.com/data/national as of July 6, 2021 -- default = 0.058204761
    R_0 = 1.66; % basic reproduction number (new infections/existing infection) from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0249271 -- default = 1.66
    value_stat_life_tot = 11.5e6/1.0; % value of a statistical life by U.S. Environmental Protection Agency and the Department of Transportation (USD) from Glover et al. (2020) -- default = 11.5e6
    c_per_capita_w_us = 45175/52; % U.S. weekly per capita consumption (USD/week = (USD/year)/(weeks/year)) from Glover et al. (2020) -- default = 45175/52
    weeks_life_remain = 37*52; % remaining life (weeks = years * weeks/year) -- default = 37*52
    % other source: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0241536

    % externally set parameters
    PAR = struct();
    PAR.INTEREST_A = 0.030; % interest rate (annual) -- default = 0.030
    PAR.ALPHA = 0.439; % Cobb-Douglas factor share of intermediate inputs -- default = 0.439
    PAR.GAMMA = 0.500; % productivity share of infected workers -- default = 0.500
    PAR.SIGMA = 1.000; % utility CRRA coefficient (1.0 = log utility) -- default = 1.000
    PAR.I_INIT = 1e2/331e6; % initial fraction of population infected -- default = 1e3/331e6
    PAR.TFP2 = 0.000; % TFP penalty from infections -- default = 0.000
    PAR.TFP3 = 0.000; % TFP curvature w.r.t. infections -- default = 0.000
    PAR.TT = 52; % weeks from start of pandemic until arrival of vaccine -- default = 52
    PAR.TTP = 104; % weeks from start of pandemic until end of simulation -- default = 104

    % grid parameters
    PAR.N_GRID_I = 200; % number of grid points for share of infected (S) -- default = 200
    PAR.N_GRID_R = 200; % number of grid points for share of infected (S) -- default = 200
    PAR.N_GRID_L = 200; % number of grid points for share of infected (S) -- default = 200
    PAR.GRID_POWER = 2; % power to which [0,1]-grids are raised -- default = 2

    % derived parameters
    PAR.INTEREST_W = (1 + PAR.INTEREST_A)^(1/52) - 1; % interest rate (weekly)
    PAR.BETA = 1/(1 + PAR.INTEREST_W); % discount factor (weekly)

    % calibration equations
    syms rho_1 rho_2 rho_3 rho_4
    eqns = [ ...
        1/(rho_3 + rho_4) == length_infection, ...
        rho_4/(rho_3 + rho_4) == mortality_rate, ...
        (rho_1 + rho_2)/(rho_3 + rho_4) == R_0, ...
        rho_1 == rho_2/2 ...
        ];
    rho_vector = solve(eqns, [rho_1 rho_2 rho_3 rho_4]);
    value_stat_life_w = value_stat_life_tot*PAR.INTEREST_W/(1 + PAR.INTEREST_W)/(1 - 1/(1 + PAR.INTEREST_W)^weeks_life_remain);

    % calibrated parameters
    PAR.TFP1 = (c_per_capita_w_us/((1 - PAR.ALPHA)*(PAR.ALPHA/PAR.INTEREST_W)^(PAR.ALPHA/(1 - PAR.ALPHA))))^(1 - PAR.ALPHA); % TFP intercept
    PAR.RHO1 = eval(rho_vector.rho_1); % SIR model parameter: power of lockdown (weekly)
    PAR.RHO2 = eval(rho_vector.rho_2); % SIR model parameter: infection rate (weekly)
    PAR.RHO3 = eval(rho_vector.rho_3); % SIR model parameter: recovery rate (weekly)
    PAR.RHO4 = eval(rho_vector.rho_4); % SIR model parameter: death rate (weekly)
    if PAR.SIGMA == 1.0
        PAR.NU = value_stat_life_w/c_per_capita_w_us - log(c_per_capita_w_us); % value of being alive
    else
        PAR.NU = value_stat_life_w*c_per_capita_w_us^(-PAR.SIGMA) - (c_per_capita_w_us^(1 - PAR.SIGMA) - 1)/(1 - PAR.SIGMA); % value of being alive
    end
    
    % parameters reset after calibration
    if l == 2
        PAR.BETA = 1 - 0.5*(1 - PAR.BETA); % low discount factor
    elseif l == 3
        PAR.BETA = 1 - 2.0*(1 - PAR.BETA); % high discount factor
    elseif l == 4
        PAR.ALPHA = 0.5*PAR.ALPHA; % low Cobb-Douglas factor share of intermediate inputs
    elseif l == 5
        PAR.ALPHA = 2.0*PAR.ALPHA; % high Cobb-Douglas factor share of intermediate inputs
    elseif l == 6
        PAR.GAMMA = 1 - 0.5*(1 - PAR.GAMMA); % low productivity penalty from being infected
    elseif l == 7
        PAR.GAMMA = 1 - 2.0*(1 - PAR.GAMMA); % high productivity penalty from being infected
    elseif l == 8
        PAR.NU = 0.5*PAR.NU; % low value of life
    elseif l == 9
        PAR.NU = 2.0*PAR.NU; % high value of life
    elseif l == 10
        PAR.RHO1 = 0.5*PAR.RHO1; % low additional transmission rate at work
    elseif l == 11
        PAR.RHO1 = 2.0*PAR.RHO1; % high additional transmission rate at work
    elseif l == 12
        PAR.RHO2 = 0.5*PAR.RHO2; % low baseline transmission rate outside of work
    elseif l == 13
        PAR.RHO2 = 2.0*PAR.RHO2; % high baseline transmission rate outside of work
    elseif l == 14
        PAR.TT = 0.5*PAR.TT; % low duration of pandemic before vaccine arrival
    elseif l == 15
        PAR.TT = 2.0*PAR.TT; % high duration of pandemic before vaccine arrival
    end


    %%% baseline simulation and binary comparative statics
    if switch_sim || (switch_comp_binary && l < 16)
        
        % baseline simulation
        [lockdown_disc, lockdown_disc_pand, lockdown_disc_year, lockdown_thresh, lockdown_thresh_pand, lockdown_thresh_year, gdp, gdp_pand, gdp_year, D, D_pand, D_year, alive_disc, alive_disc_pand, alive_disc_year, V, V_pand, V_year, cons_equiv_loss, cons_equiv_loss_pand, cons_equiv_loss_year] = simulations(PAR, DIR_SIM, [], []);

        % display results
        fprintf('\n* NPV of lockdown\n')
        disp(['   with no pandemic   = ', num2str(lockdown_disc.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_disc.commit)])
        disp(['   with no commitment = ', num2str(lockdown_disc.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_disc.no_lock)])

        fprintf('\n* NPV of lockdown during pandemic\n')
        disp(['   with no pandemic   = ', num2str(lockdown_disc_pand.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_disc_pand.commit)])
        disp(['   with no commitment = ', num2str(lockdown_disc_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_disc_pand.no_lock)])

        fprintf('\n* NPV of lockdown during first year\n')
        disp(['   with no pandemic   = ', num2str(lockdown_disc_year.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_disc_year.commit)])
        disp(['   with no commitment = ', num2str(lockdown_disc_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_disc_year.no_lock)])

        fprintf('\n* Number of weeks with lockdown above threshold level\n')
        disp(['   with no pandemic   = ', num2str(lockdown_thresh.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_thresh.commit)])
        disp(['   with no commitment = ', num2str(lockdown_thresh.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_thresh.no_lock)])

        fprintf('\n* Number of weeks with lockdown above threshold level during pandemic\n')
        disp(['   with no pandemic   = ', num2str(lockdown_thresh_pand.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_thresh_pand.commit)])
        disp(['   with no commitment = ', num2str(lockdown_thresh_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_thresh_pand.no_lock)])

        fprintf('\n* Number of weeks with lockdown above threshold level during first year\n')
        disp(['   with no pandemic   = ', num2str(lockdown_thresh_year.no_pand)])
        disp(['   with commitment    = ', num2str(lockdown_thresh_year.commit)])
        disp(['   with no commitment = ', num2str(lockdown_thresh_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(lockdown_thresh_year.no_lock)])

        fprintf('\n* NPV of GDP or wage bill or consumption\n')
        disp(['   with no pandemic   = ', num2str(gdp.no_pand/gdp.no_pand)])
        disp(['   with commitment    = ', num2str(gdp.commit/gdp.no_pand)])
        disp(['   with no commitment = ', num2str(gdp.no_commit/gdp.no_pand)])
        disp(['   with no lockdown   = ', num2str(gdp.no_lock/gdp.no_pand)])

        fprintf('\n* NPV of GDP or wage bill or consumption during pandemic\n')
        disp(['   with no pandemic   = ', num2str(gdp_pand.no_pand/gdp_pand.no_pand)])
        disp(['   with commitment    = ', num2str(gdp_pand.commit/gdp_pand.no_pand)])
        disp(['   with no commitment = ', num2str(gdp_pand.no_commit/gdp_pand.no_pand)])
        disp(['   with no lockdown   = ', num2str(gdp_pand.no_lock/gdp_pand.no_pand)])

        fprintf('\n* NPV of GDP or wage bill or consumption during first year\n')
        disp(['   with no pandemic   = ', num2str(gdp_year.no_pand/gdp_year.no_pand)])
        disp(['   with commitment    = ', num2str(gdp_year.commit/gdp_year.no_pand)])
        disp(['   with no commitment = ', num2str(gdp_year.no_commit/gdp_year.no_pand)])
        disp(['   with no lockdown   = ', num2str(gdp_year.no_lock/gdp_year.no_pand)])
        
        fprintf('\n* Accumulated share of deaths\n')
        disp(['   with no pandemic   = ', num2str(D.no_pand)])
        disp(['   with commitment    = ', num2str(D.commit)])
        disp(['   with no commitment = ', num2str(D.no_commit)])
        disp(['   with no lockdown   = ', num2str(D.no_lock)])

        fprintf('\n* Accumulated share of deaths during pandemic\n')
        disp(['   with no pandemic   = ', num2str(D_pand.no_pand)])
        disp(['   with commitment    = ', num2str(D_pand.commit)])
        disp(['   with no commitment = ', num2str(D_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(D_pand.no_lock)])

        fprintf('\n* Accumulated share of deaths during first year\n')
        disp(['   with no pandemic   = ', num2str(D_year.no_pand)])
        disp(['   with commitment    = ', num2str(D_year.commit)])
        disp(['   with no commitment = ', num2str(D_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(D_year.no_lock)])

        fprintf('\n* Discounted sum of shares alive\n')
        disp(['   with no pandemic   = ', num2str(alive_disc.no_pand)])
        disp(['   with commitment    = ', num2str(alive_disc.commit)])
        disp(['   with no commitment = ', num2str(alive_disc.no_commit)])
        disp(['   with no lockdown   = ', num2str(alive_disc.no_lock)])

        fprintf('\n* Discounted sum of shares alive during pandemic\n')
        disp(['   with no pandemic   = ', num2str(alive_disc_pand.no_pand)])
        disp(['   with commitment    = ', num2str(alive_disc_pand.commit)])
        disp(['   with no commitment = ', num2str(alive_disc_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(alive_disc_pand.no_lock)])

        fprintf('\n* Discounted sum of shares alive during first year\n')
        disp(['   with no pandemic   = ', num2str(alive_disc_year.no_pand)])
        disp(['   with commitment    = ', num2str(alive_disc_year.commit)])
        disp(['   with no commitment = ', num2str(alive_disc_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(alive_disc_year.no_lock)])

        fprintf('\n* Welfare\n')
        disp(['   with no pandemic   = ', num2str(V.no_pand)])
        disp(['   with commitment    = ', num2str(V.commit)])
        disp(['   with no commitment = ', num2str(V.no_commit)])
        disp(['   with no lockdown   = ', num2str(V.no_lock)])

        fprintf('\n* Welfare during pandemic\n')
        disp(['   with no pandemic   = ', num2str(V_pand.no_pand)])
        disp(['   with commitment    = ', num2str(V_pand.commit)])
        disp(['   with no commitment = ', num2str(V_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(V_pand.no_lock)])

        fprintf('\n* Welfare during first year\n')
        disp(['   with no pandemic   = ', num2str(V_year.no_pand)])
        disp(['   with commitment    = ', num2str(V_year.commit)])
        disp(['   with no commitment = ', num2str(V_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(V_year.no_lock)])

        fprintf('\n* Consumption-equivalent welfare loss from pandemic\n')
        disp(['   with no pandemic   = ', num2str(0)]) % cons_equiv_loss.no_pand
        disp(['   with commitment    = ', num2str(cons_equiv_loss.commit)])
        disp(['   with no commitment = ', num2str(cons_equiv_loss.no_commit)])
        disp(['   with no lockdown   = ', num2str(cons_equiv_loss.no_lock)])

        fprintf('\n* Consumption-equivalent welfare loss from pandemic during pandemic\n')
        disp(['   with no pandemic   = ', num2str(0)]) % cons_equiv_loss_pand.no_pand
        disp(['   with commitment    = ', num2str(cons_equiv_loss_pand.commit)])
        disp(['   with no commitment = ', num2str(cons_equiv_loss_pand.no_commit)])
        disp(['   with no lockdown   = ', num2str(cons_equiv_loss_pand.no_lock)])

        fprintf('\n* Consumption-equivalent welfare loss from pandemic during first year\n')
        disp(['   with no pandemic   = ', num2str(0)]) % cons_equiv_loss_year.no_pand
        disp(['   with commitment    = ', num2str(cons_equiv_loss_year.commit)])
        disp(['   with no commitment = ', num2str(cons_equiv_loss_year.no_commit)])
        disp(['   with no lockdown   = ', num2str(cons_equiv_loss_year.no_lock)])

        if PAR.SIGMA == 1.0
            fprintf('\n* Consumption-equivalent welfare loss from no commitment\n')
            disp(['   commitment vs. no commitment = ', num2str(1 - exp((V.no_commit - V.commit)/alive_disc.commit))])

            fprintf('\n* Consumption-equivalent welfare loss from no commitment during pandemic\n')
            disp(['   commitment vs. no commitment = ', num2str(1 - exp((V_pand.no_commit - V_pand.commit)/alive_disc_pand.commit))])

            fprintf('\n* Consumption-equivalent welfare loss from no commitment during first year\n')
            disp(['   commitment vs. no commitment = ', num2str(1 - exp((V_year.no_commit - V_year.commit)/alive_disc_year.commit))])
        end
    end
end


%%% multiple comparative statics
if switch_comp_multiple
    
    % define parameters and values to be looped over for comparative statics
    par_name = {'BETA', 'RHO1', 'RHO2', 'ALPHA', 'GAMMA', 'NU'};
    BETA_min = min(max(1 - (1 - PAR.BETA)*2.0, 0), 1);
    BETA_max = min(max(1 - (1 - PAR.BETA)/2, 0), 1);
    RHO1_min = PAR.RHO1/2;
    RHO1_max = PAR.RHO1*2.0;
    RHO2_min = PAR.RHO2/2;
    RHO2_max = PAR.RHO2*2.0;
    ALPHA_min = .05;
    ALPHA_max = .95;
    GAMMA_min = min(max(1 - (1 - PAR.GAMMA)*2.0, 0), 1);
    GAMMA_max = min(max(1 - (1 - PAR.GAMMA)/2, 0), 1);
    NU_min = PAR.NU/2;
    NU_max = PAR.NU*2.0;
    n_comp_stat = 2;
    par_range = { ...
        unique([linspace(BETA_min, PAR.BETA, n_comp_stat), linspace(PAR.BETA, BETA_max, n_comp_stat)]), ... % annual discount rates
        unique([linspace(RHO1_min, PAR.RHO1, n_comp_stat), linspace(PAR.RHO1, RHO1_max, n_comp_stat)]), ... % SIR model: transmission parameter at work
        unique([linspace(RHO2_min, PAR.RHO2, n_comp_stat), linspace(PAR.RHO2, RHO2_max, n_comp_stat)]), ... % SIR model: transmission parameter outside of work
        unique([linspace(ALPHA_min, PAR.ALPHA, n_comp_stat), linspace(PAR.ALPHA, ALPHA_max, n_comp_stat)]), ... % Cobb-Douglas factor share of intermediate inputs
        unique([linspace(GAMMA_min, PAR.GAMMA, n_comp_stat), linspace(PAR.GAMMA, GAMMA_max, n_comp_stat)]), ... % productivity share of infected workers
        unique([linspace(NU_min, PAR.NU, n_comp_stat), linspace(PAR.NU, NU_max, n_comp_stat)]) ... % flow value of life
    };

    % delete previous results
    FILE_COMP_STAT = [DIR_SIM, '/sim_comp_stat.xlsx'];
    if exist(FILE_COMP_STAT, 'file')
        delete(FILE_COMP_STAT)
    end
    
    % loop over parameters
    parfor vv = 1:length(par_name)
      nv = length(par_range{vv});
      lockdown_disc_no_pand = NaN(nv, 1);
      lockdown_disc_commit = NaN(nv, 1);
      lockdown_disc_no_commit = NaN(nv, 1);
      lockdown_disc_no_lock = NaN(nv, 1);
      lockdown_thresh_no_pand = NaN(nv, 1);
      lockdown_thresh_commit = NaN(nv, 1);
      lockdown_thresh_no_commit = NaN(nv, 1);
      lockdown_thresh_no_lock = NaN(nv, 1);
      gdp_no_pand = NaN(nv, 1);
      gdp_commit = NaN(nv, 1);
      gdp_no_commit = NaN(nv, 1);
      gdp_no_lock = NaN(nv, 1);
      deaths_no_pand = NaN(nv, 1);
      deaths_commit = NaN(nv, 1);
      deaths_no_commit = NaN(nv, 1);
      deaths_no_lock = NaN(nv, 1);
      alive_disc_no_pand = NaN(nv, 1);
      alive_disc_commit = NaN(nv, 1);
      alive_disc_no_commit = NaN(nv, 1);
      alive_disc_no_lock = NaN(nv, 1);
      V_no_pand = NaN(nv, 1);
      V_commit = NaN(nv, 1);
      V_no_commit = NaN(nv, 1);
      V_no_lock = NaN(nv, 1);
      is_baseline = NaN(nv, 1);

      % loop over parameter values
      for pp = 1:nv
        fprintf('Comparative statics w.r.t. %s (%10.8f)\n', par_name{vv}, par_range{vv}(pp));

        % simulation with pandemic
        [lockdown_disc, ~, ~, lockdown_thresh, ~, ~, gdp, ~, ~, deaths, ~, ~, alive_disc, ~, ~, V, ~, ~]  = simulations(PAR, DIR_SIM, par_name{vv}, par_range{vv}(pp));

        % store results 
        lockdown_disc_no_pand(pp) = lockdown_disc.no_pand;
        lockdown_disc_commit(pp) = lockdown_disc.commit;
        lockdown_disc_no_commit(pp) = lockdown_disc.no_commit;
        lockdown_disc_no_lock(pp) = lockdown_disc.no_lock;

        lockdown_thresh_no_pand(pp) = lockdown_thresh.no_pand;
        lockdown_thresh_commit(pp) = lockdown_thresh.commit;
        lockdown_thresh_no_commit(pp) = lockdown_thresh.no_commit;
        lockdown_thresh_no_lock(pp) = lockdown_thresh.no_lock;

        gdp_no_pand(pp) = gdp.no_pand;
        gdp_commit(pp) = gdp.commit;
        gdp_no_commit(pp) = gdp.no_commit;
        gdp_no_lock(pp) = gdp.no_lock;

        deaths_no_pand(pp) = deaths.no_pand;
        deaths_commit(pp) = deaths.commit;
        deaths_no_commit(pp) = deaths.no_commit;
        deaths_no_lock(pp) = deaths.no_lock;

        alive_disc_no_pand(pp) = alive_disc.no_pand;
        alive_disc_commit(pp) = alive_disc.commit;
        alive_disc_no_commit(pp) = alive_disc.no_commit;
        alive_disc_no_lock(pp) = alive_disc.no_lock;

        V_no_pand(pp) = V.no_pand;
        V_commit(pp) = V.commit;
        V_no_commit(pp) = V.no_commit;
        V_no_lock(pp) = V.no_lock;

        is_baseline(pp) = (par_range{vv}(pp) == PAR.(par_name{vv}));
      end

      % write results to Excel file
      results_comp_stat = table( ...
         par_range{vv}', ...
         lockdown_disc_no_pand, lockdown_disc_commit, lockdown_disc_no_commit, lockdown_disc_no_lock, ...
         lockdown_thresh_no_pand, lockdown_thresh_commit, lockdown_thresh_no_commit, lockdown_thresh_no_lock, ...
         gdp_no_pand, gdp_commit, gdp_no_commit, gdp_no_lock, ...
         deaths_no_pand, deaths_commit, deaths_no_commit, deaths_no_lock, ...
         alive_disc_no_pand, alive_disc_commit, alive_disc_no_commit, alive_disc_no_lock, ...
         V_no_pand, V_commit, V_no_commit, V_no_lock, ...
         is_baseline ...
      );
      results_comp_stat.Properties.VariableNames(1) = {'par_value'};
      writetable(results_comp_stat, FILE_COMP_STAT, 'Sheet', par_name{vv});
    end
end

    
%%% final housekeeping
% end parallel pool
fprintf('\n')
delete(gcp('nocreate'))

% remove parallelization packages
rmdir(DIR_PAR, 's')

% stop timer and show elapsed time
toc

% close log file
fprintf('\nDONE!\n')
diary close
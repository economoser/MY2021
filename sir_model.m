function [ss, ii, rr, dd] = sir_model(I, R, L, PAR, ind)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIR_MODEL  Return next period's health state, assuming no more new
    % infections can occur after the arrival of the vaccine. Called by
    % simulations.m, backward_induct_comm.m, backward_induct_nocomm.m,
    % gdp_vax.m, and utility_vax.m.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Inputs:
    %   I = current share of infected individuals.
    %   R = current share of recovered individuals.
    %   L = current lockdown policy.
    %   PAR = structure array of parameters.
    %   ind = true (if after arrival of the vaccine) or false (if before
    %      arrivale of the vaccine).
    %
    % Outputs:
    %   ss = share of susceptible individuals next period.
    %   ii = share of infected individuals next period.
    %   rr = share of recovered individuals next period.
    %   dd = share of dead individuals next period.
    
    
    %%% solve for state next period
    if ~ind % without vaccine
        S = 1 - I - R*(1 + PAR.RHO4/PAR.RHO3);
        ss = (1 - (PAR.RHO1*L.^2 + PAR.RHO2).*I).*S;
        ii = (1 - PAR.RHO3 - PAR.RHO4)*I + (PAR.RHO1*L.^2 + PAR.RHO2).*I.*S;
        ii(ss < 0) = ii(ss < 0) + ss(ss < 0);
        ss(ss < 0) = 0;
        rr = (1 - ss - ii)/(1 + PAR.RHO4/PAR.RHO3);
        dd = rr*PAR.RHO4/PAR.RHO3;
    else % with vaccine
        S = 1 - I - R*(1 + PAR.RHO4/PAR.RHO3);
        ss = S;
        ii = (1 - PAR.RHO3 - PAR.RHO4)*I;
        ii(ss < 0) = ii(ss < 0) + ss(ss < 0);
        ss(ss < 0) = 0;
        rr = (1 - ss - ii)/(1 + PAR.RHO4/PAR.RHO3);
        dd = rr*PAR.RHO4/PAR.RHO3;
    end
    
    
    %%% ensure that health state is feasible
    ss(S < 0) = 0; % if S < 0, then we are in an area of (I,R) that is not feasible, in which case we impose (I(t+1), R(t+1)) = (0, R(t)).
    ii(S < 0) = 0;
    rr(S < 0) = R(S < 0);
    dd(S < 0) = R(S < 0)*PAR.RHO4/PAR.RHO3;
end
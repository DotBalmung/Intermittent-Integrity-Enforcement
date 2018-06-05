function emax = PlotErrorCurveMoreF(A,B,C,D,G,select,x0,Q,R,Pfalse,ConLvl,policy,probIncrease,precision,f)

    % Plots the Error Curve - either time for specific policy or L vs emax
    % Inputs:
    %   A,B,C,D: matrices of the system
    %   G: matrix of the attackers influence, I by default
    %   select: 0-time series, 1- L vs emax
    %   x0: Starting state, 0 vector by default
    %   Q,R: Noise matrices for states and sensors, I by default
    %   Pfalse: probability of false positive of the detector, default 5%
    %   ConLvl: confidence level for variability of estimation error, ~0 by
    %       default, since we mostly care only about \delta e
    %   policy(select0): [TimeEndOfsearch howOften]
    %   policy(select1): upper bound on L (inter-enforcement distance)
    %   probIncrease: Probability of detection caused by the attacker
    %   precision: Precision for computation of alpha parameter, 10^-5 def.
    %   f: Number of encrypted points
    % Outputs:
    %   emax: vector containing error values for given time or L
    % Author: Ilija Jovanov
    %
    %   ATTENTION! - in interest of keepeing everything computationally
    %   efficient, maxSearch parameter limits number of time points where
    %   search is performed. Thus, results might get incorect for high
    %   values of inter-enforcement distance. Up to maxSearch/4 should be
    %   fine
    maxSearch = 200;
 
%% define default values for empty inputs
    n = size(A,1); %number of states
    p = size(C,1); %number of sensors
    if isempty(G)
        G = eye(p);
    end
    q = size(G,2); %number of attacked sensors
    if isempty(select)
        select = 0;
    end
    if isempty(x0)
        x0 = zeros(n,1);
    end
    if isempty(Q)
        Q = eye(n);
    end
    if isempty(R)
        R = eye(p);
    end
    if isempty(Pfalse)
        Pfalse = 0.05;
    end
    xChi = 0:0.01:100; %setup x axis for chi square. For lower order systems 100 is enough upper limit, make sure for higher order systems
    ChiCDF = chi2cdf(xChi,p); %compute cdf, to make sure what is the probability of detection for specific threshold
    h = xChi(find(ChiCDF>(1-Pfalse),1)); %find the appropriate threshold of the detector
    if isempty(ConLvl)
        ConLvl = 0.0001;
    end
    if isempty(policy)
        if select
            policy = 10; %test up to L=10
        else
            policy = [200 0]; %200 steps, no enforcment
        end
    end
    if isempty(probIncrease)
        probIncrease = 0.001;
    end
    if isempty(precision)
        precision = 5;
    end

%% Plot requested
if ~select % If time series of max(e) were requested
    locPolicy = FormGlobalPeriodicPolicyMoreF(q,min(maxSearch,policy(1)),policy(2),f);
    ErrorNorms = FindErrorRegionsCurve(A,B,C,D,G,x0,Q,R,h,ConLvl,locPolicy,probIncrease,precision);
    plot(ErrorNorms,'x-');
    xlabel('Time sample number [k]');
    ylabel('Maximum introduced error [emax]');
    emax = ErrorNorms;
else % If L vs e was requested
    emax = zeros(policy,1);
    duration = min(maxSearch,4*policy); %4x policy so attack dynamics stabilizes
    for i = 2:policy
        disp(['Computing for L = ' num2str(i)]);
        locPolicy = FormGlobalPeriodicPolicyMoreF(q,duration,i,f); %keep duration constant to keep alpha stable
        allErrors = FindErrorRegionsCurve(A,B,C,D,G,x0,Q,R,h,ConLvl,locPolicy,probIncrease,precision);
        emax(i) = max(allErrors(1:floor(duration/i)*i)); %look only at the periodic ones
    end
    plot(emax);
    xlabel('Inter-enforcement distance [L]');
    ylabel('Maximum introudced error [emax]');
end

end
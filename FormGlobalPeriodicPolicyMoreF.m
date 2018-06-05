function policy = FormGlobalPeriodicPolicyMoreF(q,duration,period,f)
    % Forms the authentication policy
    % Inputs:
    %   q: Number of compromised sensors
    %   duration: Length of the simulation
    %   period: How often the data is authenticated
    %   f: How many data points
    % Outputs:
    %   policy: Binary matrix: 1-encrypted and 0-not encrypted
    %           row # corresponds to sensor #, and column # to time k
    % Author: Ilija Jovanov
    if period == 0
        policy = zeros(q,duration);
    else
        policy = [1:duration];
        policy = mod(policy,period)==0;
        originalPolicy = policy;
        for i=2:f
            policy = policy + [originalPolicy(i:length(originalPolicy)) zeros(1,i-1)];
        end
        policy = repmat(policy,q,1);
    end
end
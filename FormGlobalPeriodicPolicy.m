function policy = FormGlobalPeriodicPolicy(q,duration,period)
    % Forms the authentication policy
    % Inputs:
    %   q: Number of compromised sensors
    %   duration: Length of the simulation
    %   period: How often the data is authenticated
    % Outputs:
    %   policy: Binary matrix: 1-encrypted and 0-not encrypted
    %           row # corresponds to sensor #, and column # to time k
    % Author: Ilija Jovanov
    if period == 0
        policy = zeros(q,duration);
    else
        policy = [1:duration];
        policy = mod(policy,period)==0;
        policy = repmat(policy,q,1);
    end
end
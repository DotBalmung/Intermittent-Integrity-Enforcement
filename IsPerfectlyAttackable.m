function attkbl=IsPerfectlyAttackable(A,C,K,Gatk)
% Description:
%	Tests wether or not the system is perfectly attackable (PA)
% Inputs:
%	A: State dynamics of the system, state-space matrix A
%   C: Outputs matrix of the system, state-space matrix C
%   K: Kalman gain in steady state
%   Gatk: Sparse column identity matrix, mapping attacks to correct places
% Outputs:
%   attkbl: -1 if not unstable (not PA), 0 if not PA, and 1 if PA
% Author: Ilija Jovanov 2017

% compute eigenvalues of A
[V,D] = eig(A);
% form a diagonal matrix of eigenvalues
lambdas = diag(D);
% find all unstable eigenvalues of A
unst_ind = find(abs(lambdas)>=1);
% If there are no unstable eigenvalues, system is not PA, and no need for
% second step (checking wether v is reachable state of dynamics from paper)
if isempty(unst_ind)
    attkbl = -1;
    %disp('system is not perfectly attackable');
else
    attkbl = 0;
    for i=1:length(unst_ind)
        % choose an eigenvector corresponding to an unstable eigenvalue
        p = V(:,unst_ind(i));
        % Compute C*v to observe it's support (refer to CDC2017 paper - PA)
        Cp = C*p;
        % Check the conditions
        % discussion on topic how they can be checked can be found at
        % https://www.mathworks.com/matlabcentral/newsreader/view_thread/270525?requestedDomain=www.mathworks.com
        Z = orth(Gatk);
        if norm(round(Cp-Z*(Z\Cp),6)) == 0 %check if Cp is in the span of Gatk. If it is, check the second condition
            Z = orth(ctrb(A-K*C*A,K*Gatk));
            if norm(round(p-Z*(Z\p),6)) == 0 %check if p belongs to the span of the controllability matrix of (A-KCA,KGatk)
                attkbl = 1;
            end
        end
    end
    % Write in console whether or not the system is attackable
    if attkbl == 0
        %disp('system is not perfectly attackable');
    else
        %disp('system is perfectly attackable');
    end
end

end
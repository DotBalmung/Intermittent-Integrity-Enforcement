function [ErrorNorms,attVecs] = PlanAttackGAasEA_CDC18(A,C,G,K,iterations,encT)
% Description:
%	Plans the attack based on CDC2017 paper method. Modified to provide
%	quick results, although suboptimal. Essentially, we assume that the
%	attacker attacks using greedy attack, and then returns along the same
%	path. Since this can sometimes cross stealthy requirements, the attack
%	is scaled appropriately to compensate, and keep the attack stealthy.
% Inputs:
%	A: State-space matrix of the system (A)
%   C: State-space matrix of the system (C)
%   G: Matrix that maps attacks into corresponding sensors
%   K: Kalman gain
%   iterations: number of iterations for which attack is planned
%   encT: -1 indicates attack planning with no authentication and looks
%   "iterations" steps ahead, 0 utilizes greedy attack (looks only one step
%   ahead), and encT>0 looks encT steps ahead when planning attack
% Outputs:
%	ErrorNorms: Vector containing all 2-norms of estimation error for all k
%   attVecs: matrix containing attack vectors used to reach this error
% Author: Ilija Jovanov 2016

%[~,a2] = PlanAttack(A,C,G,K,floor(encT/2),0); % greedy attack planning
[~,a2] = PlanAttack(A,C,G,K,encT,0); % greedy attack planning
%a2 = [a2 1*a2(:,ceil(encT/2)-1:-1:1) zeros(size(a2,1),1)]; % form a greedy attack that returns to 0 along the same path
%a2 = [a2 zeros(size(a2,1),1)]; % form a greedy attack that returns to 0 along the same path

%% design attack
times = ceil(iterations/encT);
a3 = repmat(a2,1,2*times); % repeat increase-decrease of the attack as required to keep the attack on all data
a3 = a3/max(max(abs(ComputeZ(A,C,K,G,a3)))); % scale with highest obtained residue, so that attack remains stealthy

%[z3,e3] = ComputeZ(A,C,K,G,a3); % used for testing of the function during development

a3 = a3(:,1:iterations); % truncate to required length
%z3 = z3(:,1:iterations);
%e3 = e3(:,1:iterations);

%ErrorNorms = e3; % Originally used for testing purposes, so I've left it in code, so I can test potential future modifications
ErrorNorms = 0; % Not really used anymore, so I've left it as 0. This way it's easier to notice if I attempt to use e3 in codes by mistake.
attVecs = a3;
%residues = z3;

end
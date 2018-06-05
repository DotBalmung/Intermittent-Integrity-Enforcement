function [ErrorNorms] = FindErrorRegionsCurve(A,B,C,D,G,x0,Q,R,h,ConLvl,policy,probIncrease,precision)
  %Barebones version - does not return whole region, just curves, as it's
  %more computationally efficient this way
    %This function accepts the information about the system and the
    %integrity enforcement policy, and provides the information about the
    %regions of the attack. This function is written for DISCRETE TIME
    %SYSTEMS. As such, please discretize your system before using this
    %function.
    %Inputs:
    %   A,B,C,D: Matrices that describe the state-space of the system
    %   G: Matrix that maps the attacks into proper sensor positions
    %   x0: Starting state of the system
    %   Q,R: Parameters describing the LQ optimal controller
    %   h: Detection level for the detector function (gk)
    %   ConLvl: Confidence level - controls the probability that e_k^a is
    %           inside of the obtained regions
    %   policy: Binary matrix that describes which (compromised) sensors are 
    %       encrypted and when. Rows correspond to sensors, and columns 
    %       correspond to time steps.
    %   probIncrease: increase in the probability of he detection
    %       introduced by the attacker.
    %   precision: how precisely the threshold will be computed. Affects
    %       the precision of the obtained regions/limits. Precision will be
    %       up to 10^(-precision).
    %Outputs:
    %   ErrorNorms: Array of maximal error norms for each region
    %Author: Ilija Jovanov
    
%% Preliminary computations
% disp('Running FindErrorRegions function...');
% disp('Performing Preliminary computations...');
% Find kalman gain, Error covariance matrix, and residue covariance matrix
[K,P,Pk] = findK(A,C,Q,R,x0);
% Display whether the system is perfectly attackable
% if IsPerfectlyAttackable(A,C,K,G)
%     disp('This system is Perfectly attackable');
% else
%     disp('This system is NOT Perfectly attackable');
% end
% added as in paper, so that it's not simulation dependant
    P = C*Pk*C'+eye(size(C,1)); %unit noise matrix
% Define parameters required for the future calculations
Sig = inv(P);  % Used later in the code in several places, saving computation time in this way
p = size(C,1); % Number of sensors
q = size(G,2); % Number of compromised sensors
n = size(A,1); % Number of states of the system
CovCoeff = ConLvl*Pk;
steps = size(policy,2); % Duration of computation. Depends on the information provided about the policy

%% Compute conditions for each k-reachable region
% disp('Computing limits for each of the regions...');
% Find initial matrices for the method
CumPol = policy(:,1); %Keeps the track of all a_k^i=0 over time
M{1} = K*G; %Matrix corresponding to estimation error 
Mred{1} = M{1}(:,~CumPol); %Reduced matrix M, used in region computation
N{1} = G; %Matrix corresponding to the residue
% Find remaining conditions
for i=2:steps
    CumPol = [CumPol;policy(:,i)];
    M{i} = [(A-K*C*A)^(i-1)*K*G M{i-1}];
    Mred{i} = M{i}(:,~CumPol);
    N{i} = [-C*A*M{i-1} G];
end
% Find alpha threshold defined in the arxiv paper
alph = FindAlpha(probIncrease,steps*p,2*h+steps*p,precision);
% Form matrix theta as in the arxiv paper
The = [N{1} zeros(p,q*(steps-1))]'*Sig*[N{1} zeros(p,q*(steps-1))];
for i = 2:steps
    The = The + [N{i} zeros(p,q*(steps-i))]'*Sig*[N{i} zeros(p,q*(steps-i))];
end
% Perform reduction of theta matrix
TheRed = The(~policy,~policy);
for i = 1:steps
    RegLimits{i} = alph^2*[Mred{i} zeros(n,size(TheRed,1)-size(Mred{i},2))]/TheRed*[Mred{i} zeros(n,size(TheRed,1)-size(Mred{i},2))]' + CovCoeff;
end

%% Find the regions and error norms - this needs a better solution
% disp('Finding the regions...');
ErrorNorms = zeros(steps,1);
for k=1:steps
    ErrorNorms(k) = 1./sqrt(min(eig(inv(RegLimits{k}))));
end

end
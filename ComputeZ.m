function [z,e] = ComputeZ(A,C,K,G,a)
% Description:
%	Computes residues caused by the given attack
% Inputs:
%	A: State-space matrix of the system (A - state dynamics)
%   B: State-space matrix of the system (C - sensor dynamics)
%   K: Kalman gain
%   G: Matrix that maps attacks into corresponding sensors
%   a: attack vector
% Outputs:
%	z: residue of the system over time
%   e: state estimation error vector over time
% Author: Ilija Jovanov 2016

    k = length(a); % duration of the attack
    n = size(A,1); % number of states
    p = size(C,1); % number of sensors
    q = size(G,2); % number of compromized sensors
    
    % form non-recursive computational matrices as in CDC2017 paper
    M{1} = K*G;
    N{1} = G;
    AllNs = N{1};
    AllMs = M{1};
    % Find z and e for every attack vector
    for i=2:k
        M{i} = [(A-K*C*A)*M{i-1} K*G];
        N{i} = [-C*A*M{i-1} G];
        AllNs = [AllNs zeros((i-1)*p,q);N{i}];
        AllMs = [AllMs zeros((i-1)*n,q);M{i}];
    end
    % reshape such that columns indicate time instance k
    z = reshape(AllNs*reshape(a,[],1),p,[]);
    e = reshape(AllMs*reshape(a,[],1),n,[]);
end
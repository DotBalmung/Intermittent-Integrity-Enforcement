function [K,P,Pk] = findK(A,C,Q,R,x0)
% Description:
%	Computes Kalman gain, and appropriate covariance matrices
% Inputs:
%	A: State dynamics of the system, state-space matrix A
%   C: Outputs matrix of the system, state-space matrix C
%   Q: State noise covariance matrix of the system
%   R: Sensor noise covariance matrix of the system
%   x0: Initial system states
% Outputs:
%	K: Kalman gain
%   P: Covariance of the residue
%   Pk: Estimation error covariance matrix
% Author: Ilija Jovanov 2017

    Pk=eye(length(A)); %Sigma in paper
    actualXm1 = x0; %real state of the system (without the noise)
    xm1 = x0;
    Pkm1 = Pk;
    %iterate in order for Kalman gain K to stabilize
    for i = 1:10000
        %calculation of actual state of the system and sensor values
        actualX = A * actualXm1 + Q*randn(length(A),1);%calculate current state from last state
        actualY = C * actualX; %calculate value that sensor returns
        actualXm1 = actualX; %save current X as X(t-1) for next iteration
        %kalman filtering
        yk = actualY + R*randn(size(C,1),1);%noisy sensor measurement in this iteration
        xm = A * xm1; %state from model
        Pkm = A * Pkm1 * A' + Q; % state prediction covariance
        K = Pkm * C' / (C * Pkm * C' + R);
        xk = xm + K * (yk - C * xm); %approx. state with respect to mesurements
        z(i,:) = yk - C * xm; % compute residue
        Pk = (eye(length(A)) - K * C) * Pkm; %updated state covariance
        %update for next iteration
        xm1 = xk;
        Pkm1 = Pk;
    end
    P = cov(z);
end
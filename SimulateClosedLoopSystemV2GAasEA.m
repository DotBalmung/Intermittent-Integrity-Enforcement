function [x,y,ya,x_hat,u,a] = SimulateClosedLoopSystemV2GAasEA_CDC18(A,B,C,Gatk,W,U,Q,R,x0,IsAttk,AtkDuration,AttkScale,reference,encT)
% Desc: Simulates the system in closed loop
% Inputs:
%   A,B,C: matrices describing dynamics of the system
%   Gatk: Matrix mapping attack vector into sensors
%   W,U: Parameters for LQR state and inputs respectively
%   Q,R: Variance matrices for the state and sensor noise
%   x0: Initial state of the system
%   IsAttk: 1-there is attack, 0-no attack
%   AtkDuration: number of steps for the attack
%   AttkScale: scaling factor for the attack
%   reference: Reference signal for the controller to follow
%   encT: How often the data integrity is being enforced
% Outputs:
%   x: actual states of the system
%   y: actual sensor values of the system
%   ya: recorded sensor values of the system
%   x_hat: estimated states of the system
%   u: control signal from the LQR controller
%   a: attack vector of the system
% Author: Ilija Jovanov 2017

K = findK(A,C,Q,R,x0);
if (IsPerfectlyAttackable(A,C,K,Gatk))
    disp('System is perfectly attackable');
else
    disp('System is not perfectly attackable');
end
[~,a] = PlanAttackGAasEA(A,C,Gatk,K,AtkDuration,encT);
a = AttkScale*[zeros(size(a)) a]; %half the time no attack, half the time attack
n = size(A,1); % #states
m = size(B,2); % #inputs
p = size(C,1); % #outputs

% Initialize system parameters for the simulation
x(1:n,1) = x0; % system states
y(1:p,1) = C*x(:,1) + R*randn(p,1); % original sensor values
ya(1:p,1) = y(:,1) + IsAttk*Gatk*a(:,1); % sensor values under attack
x_hat(1:n,1) = x(1:n,1); % estimator output
u(1:m,1) = 0; % control signal (input)
[S,~,~] = dare(A,B,W,U); 
L = -(B'*S*B + U)\B'*S*A; % Gain of the controller

% Iterate to simulate obtained values over time
for i=2:AtkDuration*2
    x(:,i) = A*x(:,i-1) + B*u(:,i-1) + Q*randn(n,1); % update model states
    y(:,i) = C*x(:,i) + R*randn(p,1); % update model sensors
    ya(:,i) = y(:,i) + IsAttk*Gatk*a(:,i); % update sensors under attack
    x_hat(:,i) = A*x_hat(:,i-1)+B*u(:,i-1) + K*(ya(:,i) - C*(A*x_hat(:,i-1)+B*u(:,i-1))); % compute state estimation for steady-state kalman filter
    u(:,i) = L*(x_hat(:,i)-reference(:,i)); % compute control inputs for next iteration
end
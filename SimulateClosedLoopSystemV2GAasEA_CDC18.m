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
[~,a] = PlanAttackGAasEA_CDC18(A,C,Gatk,K,AtkDuration,encT);
a = AttkScale*[zeros(size(a)) a]; %half the time no attack, half the time attack
n = size(A,1); % #states
m = size(B,2); % #inputs
p = size(C,1); % #outputs



% Iterate to simulate obtained values over time
for j=1:AtkDuration*2/encT
    for i=j*encT-encT+1:j*encT
        if (j+i==2)
            % Initialize system parameters for the simulation
            x(1:n,1) = x0; % system states
            y(1:p,1) = C*x(:,1) + R*randn(p,1); % original sensor values
            ya(1:p,1) = y(:,1) + IsAttk*Gatk*a(:,1); % sensor values under attack
            x_hat(1:n,1) = x(1:n,1); % estimator output
            u(1:m,1) = 0; % control signal (input)
            %
            x_noAttk(1:n,1) = x0; % system states
            y_noAttk(1:p,1) = C*x(:,1) + R*randn(p,1); % original sensor values
            ya_noAttk(1:p,1) = y(:,1) + IsAttk*Gatk*a(:,1); % sensor values under attack
            x_hat_noAttk(1:n,1) = x(1:n,1); % estimator output
            u_noAttk(1:m,1) = 0; % control signal (input)
            %
            [S,~,~] = dare(A,B,W,U); 
            L = -(B'*S*B + U)\B'*S*A; % Gain of the controller
        else
            x_noAttk(:,i) = A*x_noAttk(:,i-1) + B*u_noAttk(:,i-1) + Q*randn(n,1); % update model states
            y_noAttk(:,i) = C*x_noAttk(:,i) + R*randn(p,1); % update model sensors
            ya_noAttk(:,i) = y_noAttk(:,i); % update sensors under attack
            x_hat_noAttk(:,i) = A*x_hat_noAttk(:,i-1)+B*u_noAttk(:,i-1) + K*(ya_noAttk(:,i) - C*(A*x_hat_noAttk(:,i-1)+B*u_noAttk(:,i-1))); % compute state estimation for steady-state kalman filter
            u_noAttk(:,i) = L*(x_hat_noAttk(:,i)-reference(:,i)); % compute control inputs for next iteration
            %
            x(:,i) = A*x(:,i-1) + B*u(:,i-1) + Q*randn(n,1); % update model states
            y(:,i) = C*x(:,i) + R*randn(p,1); % update model sensors
            ya(:,i) = y(:,i) + IsAttk*Gatk*a(:,i); % update sensors under attack
            x_hat(:,i) = A*x_hat(:,i-1)+B*u(:,i-1) + K*(ya(:,i) - C*(A*x_hat(:,i-1)+B*u(:,i-1))); % compute state estimation for steady-state kalman filter
            u(:,i) = L*(x_hat(:,i)-reference(:,i)); % compute control inputs for next iteration
            if i==j*encT
                x(:,i) = x_noAttk(:,i);
                y(:,i) = y_noAttk(:,i);
                ya(:,i) = ya_noAttk(:,i);
                x_hat(:,i) = x_hat_noAttk(:,i);
                u(:,i) = u_noAttk(:,i);
            end
        end
    end
end
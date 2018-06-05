function [ErrorNorms,attVecs] = PlanAttack(A,C,Gatk,K,iterations,encT)
% Description:
%	Plans the attack based on CDC2017 paper method. Forms matrices as
%	described in the paper, and then performs optimization.
% Inputs:
%	A: State-space matrix of the system (A)
%   C: State-space matrix of the system (C)
%   Gatk: Matrix that maps attacks into corresponding sensors
%   K: Kalman gain
%   iterations: number of iterations for which attack is planned
%   encT: -1 indicates attack planning with no authentication and looks
%   "iterations" steps ahead, 0 utilizes greedy attack (looks only one step
%   ahead), and encT>0 looks encT steps ahead when planning attack
% Outputs:
%	ErrorNorms: Vector containing all 2-norms of estimation error for all k
%   attVecs: matrix containing attack vectors used to reach this error
% Author: Ilija Jovanov 2016
% 07/16/2016 added a return value for the attack vectors

%% initialize matrices for calculation
    M = -(A-K*C*A);
    N = -C*A;
    P = K*Gatk;
    Q = -Gatk;
    n = size(A,1);
    m = size(C,1);
    bad = size(Gatk,2);
    ErrorNorms = zeros(iterations,1);
    attVecs = zeros(bad,iterations);
    if encT>0
        attVecs = zeros(bad,iterations*encT);
    end
    
%% form optimization problem BIG*vars=b;    
    if (encT==0)
        %no authentication case
        %form matrices
        BIGmain = [eye(n) zeros(n,n+m+bad);...
                M eye(n) zeros(n,m) P;...
                N zeros(m,n) eye(m) Q];
        BIGsupp = [zeros(m,2*n) eye(m) zeros(m,bad);...
                zeros(m,2*n) -eye(m) zeros(m,bad)];
        e0 = zeros(n,1);   
        %calculate iteratively
        for i=1:iterations
            b = [e0;zeros(n+m,1);-ones(size(BIGmain,2)-size(BIGmain,1),1)];
            %choose all possible combinations that make BIG square matrix and
            %find maximum of all of these
            combs = nchoosek(1:2*m,size(BIGmain,2)-size(BIGmain,1));
            Emax = e0;
            EnormMax = 0;
            yam = zeros(bad,1);
            for j = 1:length(combs)
                BIG = [BIGmain;BIGsupp(combs(j,:),:)];
                %check if it has full rank and if all z-s are satisfied
                if rank(BIG) == size(BIGmain,2)
                    sol = BIG\b;
                    Zvals = sol(2*n+1:2*n+m);
                    if (sum(abs(Zvals)>1)==0)
                        e1 = sol(n+1:2*n);
                        if (norm(e1)>EnormMax)
                            Emax = e1;
                            EnormMax = norm(e1);
                            yam = sol(2*n+m+1:2*n+m+bad);
                            attVecs(:,i) = sol(length(sol)-bad+1:length(sol));
                        end
                    end
                end
            end
            ErrorNorms(i,1) = EnormMax;
            e0 = Emax;
            %yam
            norm(sol(2*n+1:2*n+m));
        end
        
    elseif encT==-1
        
        activeConst = zeros(iterations*m,1);
        T = iterations;
        BIGmain = [eye(n) zeros(n,T*n+T*m+T*bad);...
                [kron(eye(T),M) zeros(T*n,n)]+[zeros(T*n,n) kron(eye(T),eye(n))] zeros(T*n,T*m) kron(eye(T),P);...
                kron(eye(T),N) zeros(T*m,n) eye(T*m) kron(eye(T),Q)];
        BIGsupp = [zeros(T*m,n+T*n) eye(T*m) zeros(T*m,T*bad);...
                zeros(T*m,n+T*n) -eye(T*m) zeros(T*m,T*bad)];
        e0 = zeros(n,1);
        b = [e0;zeros(T*n+T*m,1);-ones(size(BIGmain,2)-size(BIGmain,1),1)];
        %choose all possible combinations that make BIG square matrix and
        %find maximum of all of these
        combs = nchoosek(1:2*T*m,size(BIGmain,2)-size(BIGmain,1));
        Emax = e0;
        EnormMax = 0;
        yam = zeros(bad,1);
        AllENorms = zeros(T,1);
        for j = 1:length(combs)
            BIG = [BIGmain;BIGsupp(combs(j,:),:)];
            %check if it has full rank and if all z-s are satisfied
            if rank(BIG) == size(BIGmain,2)
                sol = BIG\b;
                Zvals = sol(n+T*n+1:n+T*n+T*m);
                if (sum(abs(Zvals)>1)==0)
                    eT = sol(T*n+1:n+T*n);
                    if (norm(eT)>EnormMax)
                        Emax = eT;
                        EnormMax = norm(eT);
                        for k = 1:T
                            AllENorms(k,1) = norm(sol(k*n+1:k*n+n));
                        end
                        activeConst = sol(T*n+n+1:T*n+n+T*m);
                        yam = sol(n+T*n+T*m+T*bad-bad+1:n+T*n+T*m+T*bad);
                    end
                end
            end
        end
        activeConst
        ErrorNorms = AllENorms;
            
    else
        %authenticate every T-th step
        %form matrices
        T = encT;
        BIGmain = [eye(n) zeros(n,T*n+T*m+T*bad);...
                [kron(eye(T),M) zeros(T*n,n)]+[zeros(T*n,n) kron(eye(T),eye(n))] zeros(T*n,T*m) kron(eye(T),P);...
                kron(eye(T),N) zeros(T*m,n) eye(T*m) kron(eye(T),Q);...
                zeros(bad,n+T*n+T*m+T*bad-bad) eye(bad)];
        BIGsupp = [zeros(T*m,n+T*n) eye(T*m) zeros(T*m,T*bad);...
                zeros(T*m,n+T*n) -eye(T*m) zeros(T*m,T*bad)];
        e0 = zeros(n,1);   
        for i=1:iterations
            disp(['Attack planning iteration ' num2str(i)]);
            b = [e0;zeros(T*n+T*m,1);zeros(bad,1);-ones(size(BIGmain,2)-size(BIGmain,1),1)]; %special zeros for last ya
            % vector is [e0,e1,...,eT,z1,...,zT,y1,...yT]^T
            %choose all possible combinations that make BIG square matrix and
            %find maximum of all of these
            combs = nchoosek(1:2*T*m,size(BIGmain,2)-size(BIGmain,1));
            Emax = e0;
            EnormMax = 0;
            yam = zeros(bad,1);
            AllENorms = zeros(T,1);
            for j = 1:length(combs)
                if ~mod(j,(length(combs)-mod(length(combs),10))/10)
                    disp(['executing... ' num2str(ceil(j/length(combs)*100)) '%']);
                end
                BIG = [BIGmain;BIGsupp(combs(j,:),:)];
                %check if it has full rank and if all z-s are satisfied
                if rank(BIG) == size(BIGmain,2)
                    sol = BIG\b;
                    Zvals = sol(n+T*n+1:n+T*n+T*m);
                    if (sum(abs(Zvals)>1)==0)
                        %eT = sol(T*n+1:n+T*n);
                        eT = sol(1:n+T*n);
                        if (norm(eT)>EnormMax)
                            Emax = eT;
                            EnormMax = norm(eT);
                            for k = 1:T
                                AllENorms(k,1) = norm(sol(k*n+1:k*n+n));
                            end
                            yam = sol(n+T*n+T*m+T*bad-bad+1:n+T*n+T*m+T*bad);
                            for k = 1:T
                                attVecs(:,(i-1)*T+k) = sol(length(sol)-(T-k+1)*bad+1:length(sol)-(T-k)*bad);
                            end
                        end
                        %MODDED - used for some testing, not in papers,
                        %allowed us to look at individual error norms
%                         for k = 1:T
%                         if (norm(sol(k*n+1:k*n+n))>EnormMax)
%                             Emax = eT;
%                             EnormMax = norm(sol(k*n+1:k*n+n));
%                             for kk = 1:T
%                                 AllENorms(kk,1) = norm(sol(kk*n+1:kk*n+n));
%                             end
%                             yam = sol(n+T*n+T*m+T*bad-bad+1:n+T*n+T*m+T*bad);
%                             for kk = 1:T
%                                 attVecs(:,(i-1)*T+kk) = sol(length(sol)-(T-kk+1)*bad+1:length(sol)-(T-kk)*bad);
%                             end
%                         end
%                         end
                        %EOMOD
                    end
                end
            end
            ErrorNorms(T*(i-1)+1:T*i,1) = AllENorms;
            e0 = Emax(T*n+1:T*n+n);
            %yam
        end
    end

end
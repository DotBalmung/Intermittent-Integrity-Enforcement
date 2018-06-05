function [ErrorNorms,attVecs] = PlanAttackGAasEAchangingL(ListOfEncs,CruiseControl)
% Description:
%	Plans attacks similarly to PlanAttackGAasEA. This function has additional
%	capability of supporting changing L values by using predefined attack
%	segments (obtained via PlanAttackGAasEA.m).
% Inputs:
%	ListOfEncs: List of inter-authentications dist.(L) provided as a vector
%   CruiseControl: Selector - 1 means cruise control case study, while in
%   the case 0 this function uses attacks for lane position tracking
% Outputs:
%	ErrorNorms: Vector containing all 2-norms of estimation error for all k
%   attVecs: matrix containing attack vectors used to reach this error
% Author: Ilija Jovanov 2017

    ErrorNorms = 0;
    
    if CruiseControl
        a = zeros(3,0);
    else
        a = zeros(4,0);
    end

    % Load all attack patterns to be used in attack planning
    CCaL1 = load('CCaLis1.mat');
    CCa{1} = CCaL1.aLis1;
    CCaL2 = load('CCaLis2.mat');
    CCa{2} = CCaL2.aLis2;
    CCaL3 = load('CCaLis3.mat');
    CCa{3} = CCaL3.aLis3;
    CCaL4 = load('CCaLis4.mat');
    CCa{4} = CCaL4.aLis4;
    CCaL5 = load('CCaLis5.mat');
    CCa{5} = CCaL5.aLis5;
    CCaL6 = load('CCaLis6.mat');
    CCa{6} = CCaL6.aLis6;
    
    LPaL1 = load('LPaLis1.mat');
    LPa{1} = LPaL1.aLis1;
    LPaL2 = load('LPaLis2.mat');
    LPa{2} = LPaL2.aLis2;
    LPaL3 = load('LPaLis3.mat');
    LPa{3} = LPaL3.aLis3;
    LPaL4 = load('LPaLis4.mat');
    LPa{4} = LPaL4.aLis4;
    
    % Form the attack sequence using attack patterns
    for i = 1:length(ListOfEncs)
        if CruiseControl
           a = [a CCa{ListOfEncs(i)}];
        else
           a = [a LPa{ListOfEncs(i)}];
        end
    end
    
    %size(a)
    attVecs = repmat(a,1,10);
    %size(attVecs)

end
function alpha=FindAlpha(e,k,h,prec)
% Description:
%	Finds the threshold for given chi^2 function and epsilon
% Inputs:
%	e: epsilon threshold
%   k: order of the chi^2 distribution
%   h: threshold of g_k function
%   prec: digits of precision
% Outputs:
%	alpha: new threshold for detection
% Author: Ilija Jovanov 2017

    beta = 1-chi2cdf(h,k);
    searchLim = 10^(-prec);
    %
    lam = 0;
    sInt = 0;
    eInt = 1000;
    % Perform bisection method to find coresponding alpha
    while searchLim<eInt-sInt
        if ((e-(1-ncx2cdf(h,k,sInt))+beta)*(e-(1-ncx2cdf(h,k,eInt/2))+beta)<0)
            sInt = sInt;
            eInt = (sInt+eInt)/2;
        else
            sInt = (sInt+eInt)/2;
            eInt = eInt;
        end
    end
    alpha = sqrt((sInt+eInt)/2);
end
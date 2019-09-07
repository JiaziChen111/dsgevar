function [lik, SIGMAu, PHI] = loglikelihoodDSGEVAR(lambda,t,p,n,GXX,GYY,GYX, mXX,mYY, mYX)
k=p*n+1;

%Compute minus log likelihood according to (A.2), DS (2004) when the dsge prior weight is finite  
if ~isinf(lambda)% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is finite.
    tmp0 = lambda*t*GYY + mYY;
    tmp1 = lambda*t*GYX + mYX;
    tmp2 = inv(lambda*t*GXX+mXX);
    SIGMAu = tmp0 - tmp1*tmp2*tmp1'; 
    SIGMAu = SIGMAu / (t*(1+lambda));
    PHI = tmp2*tmp1'; 
    prodlng1 = sum(gammaln(.5*((1+lambda)*t-n*p +1-(1:n)')));
    prodlng2 = sum(gammaln(.5*(lambda*t-n*p+1-(1:n)')));
    lik = .5*n*log(det(lambda*t*GXX+mXX))+ .5*((lambda+1)*t-k)*log(det((lambda+1)*t*SIGMAu)) ...
          - .5*n*log(det(lambda*t*GXX)) ...
          - .5*(lambda*t-k)*log(det(lambda*t*(GYY-GYX*inv(GXX)*GYX'))) ...
          + .5*n*t*log(2*pi)  ...
          - .5*log(2)*n*((lambda+1)*t-k) ...
          + .5*log(2)*n*(lambda*t-k) ...
          - prodlng1 + prodlng2;
else% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is infinite.
    iGXX = inv(GXX);
    SIGMAu = GYY - GYX*iGXX*transpose(GYX);
    PHI = iGXX*transpose(GYX);
    lik = t*(log(det(SIGMAu)) + n*log(2*pi) +  ...
                   trace(inv(SIGMAu)*(mYY - transpose(mYX*PHI) - mYX*PHI + transpose(PHI)*mXX*PHI)/t));
    lik = .5*lik;% Minus likelihood
end
end


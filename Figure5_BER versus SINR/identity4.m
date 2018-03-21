function [ p1,p2,p3,p4 ] = identity4( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e,aL)
%   lamda1 = 0.5;
%   lamda2 = 0.5;
  p1 = 0.5;
  p2 = 0.5;
  p3 = 0.5; 
  p4=0.5;
  pS = p1;
  pR = p2;
  pO = p3;
  pQ = p4;
  tau1 = 0.0001;
  tau2 = 0.0001;
  tau3 = 0.0001;
  tau4 = 0.0001;
  lamda1 = aL;
  lamda2 = aL;
  lamda3 = aL;
  lamda4 = aL;
   g1 = ( 1 - beta1 ) * trace( H11 * H11') * pS - e;
   g2 = ( 1 - beta2 ) * trace( H22 * H22') * pR - e;
   g3 = ( 1 - beta3 ) * trace( H33 * H33') * pO - e;
   g4 = ( 1 - beta4 ) * trace( H44 * H44') * pQ - e;
   lamdaS = lamda1 - tau1 * g1;
   lamdaR = lamda2 - tau2 * g2;
   lamdaO = lamda3 - tau3 * g3;
   lamdaQ = lamda4 - tau4 * g4;
       if lamdaS < 0
            lamdaS1 = 0;
         else
            lamdaS1 = lamdaS;
       end
       if lamdaR< 0
            lamdaR1 = 0;
         else
            lamdaR1 = lamdaR;
       end 
       if lamdaO< 0
            lamdaO1 = 0;
         else
            lamdaO1 = lamdaO;
       end 
        if lamdaQ< 0
            lamdaQ1 = 0;
         else
            lamdaQ1 = lamdaQ;
       end 
  lamda1 = lamdaS1;
  lamda2 = lamdaR1;
  lamda3 = lamdaO1;
  lamda4 = lamdaQ1;
 
  [ muS,muR,muO,muQ ] = get_mu4( p1,p2,p3,p4,sigma_n,beta1,beta2,beta3,beta4,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,P1_max,P2_max,P3_max,P4_max);
   gammaS = lamda1 * ( 1 - beta1 ) *  trace( H11 * H11');
   gammaR = lamda2 * ( 1 - beta2 ) *  trace( H22 * H22');
   gammaO = lamda3 * ( 1 - beta3 ) *  trace( H33 * H33');
   gammaQ = lamda4 * ( 1 - beta4 ) *  trace( H44 * H44');
  R1=( muS + gammaS )^(-1) - ( beta1 * trace( H21 * H21') * pR + beta1 * trace( H31 * H31') * pO+ beta1 * trace( H41 * H41') * pQ   + sigma_n ) / ( beta1 * trace( H11 * H11') );
     if R1 <= 0
       pS = 0;
    else 
       if R1 >= P1_max
          pS = P1_max;
       else 
          pS = R1;
       end
     end
    p1 = pS;
    
      R2=( muR + gammaR )^(-1) - ( beta1 * trace( H12 * H12') * pS + beta1 * trace( H32 * H32') * pO + beta1 * trace( H42 * H42') * pQ  + sigma_n ) / ( beta1 * trace( H22 * H22') );
     if R2 <= 0
       pR = 0;
    else 
       if R2 >= P2_max
          pR = P2_max;
       else 
          pR = R2;
       end
     end
    p2 = pR;
    
      R3=( muO + gammaO )^(-1) - ( beta1 * trace( H13 * H13') * pS + beta1 * trace( H23 * H23') * pR +  beta1 * trace( H43 * H43') * pQ  + sigma_n ) / ( beta1 * trace( H33 * H33') );
     if R3 <= 0
       pO = 0;
    else 
       if R3 >= P3_max
          pO = P3_max;
       else 
          pO = R3;
       end
     end
    p3 = pO;
      R4=( muQ + gammaQ )^(-1) - ( beta1 * trace( H14 * H14') * pS + beta1 * trace( H24 * H24') * pR + beta1 * trace( H34 * H34') * pO  + sigma_n ) / ( beta1 * trace( H44 * H44') );
     if R4 <= 0
       pQ = 0;
    else 
       if R4 >= P4_max
          pQ = P4_max;
       else 
          pQ = R4;
       end
     end
    p4 = pQ;

    fprintf( 'p1 = %3.2f; p2 = %3.2f;...\n',p1,p2,p3,p4);

     
end


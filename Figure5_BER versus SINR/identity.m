function [ p1,p2,p3] = identity(H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,P1_max,P2_max,P3_max,e,aL)
  p1 = 0.5;
  p2 = 0.5;
  p3 = 0.5; 

  pS = p1;
  pR = p2;
  pO = p3;

  tau1 = 0.01;
  tau2 = 0.01;
  tau3 = 0.01;
 
  lamda1 = aL;
  lamda2 = aL;
  lamda3 = aL;
  
  
 
    [ muS,muR,muO ] = get_mu3(  p1,p2,p3,sigma_n,beta1,beta2,beta3,H11,H22,H33,H12,H13,H21,H23,H31,H32,P1_max,P2_max,P3_max);
   gammaS = lamda1 * ( 1 - beta1 ) *  trace( H11 * H11');
   gammaR = lamda2 * ( 1 - beta2 ) *  trace( H22 * H22');
   gammaO = lamda3 * ( 1 - beta3 ) *  trace( H33 * H33');
   
  R1=( muS + gammaS )^(-1) - ( beta1 * trace( H21 * H21') * pR+ beta1 * trace( H31 * H31') * pO+ sigma_n ) / ( beta1 * trace( H11 * H11') );
     if R1 <= 0
       pS = 0;
    else 
       if R1 >= P1_max
          pS = P1_max;
       else 
          pS = R1;
       end
     end
   
    
      R2=( muR+ gammaR )^(-1) - ( beta2 * trace( H12 * H12') * pS + beta2 * trace( H32 * H32') * pO +  sigma_n ) / ( beta2 * trace( H22 * H22') );
     if R2 <= 0
       pR = 0;
    else 
       if R2 >= P2_max
          pR = P2_max;
       else 
          pR = R2;
       end
     end
   
    
      R3=( muO + gammaO)^(-1) - ( beta3 * trace( H13 * H13') * pS + beta3 * trace( H23 * H23') * pR  + sigma_n ) / ( beta3 * trace( H33 * H33') );
     if R3 <= 0
       pO = 0;
    else 
       if R3 >= P3_max
          pO = P3_max;
       else 
          pO = R3;
       end
     end
      p1 = pS;
      p2 = pR;
      p3 = pO;
      
       fprintf( 'p1 = %3.2f; p2 = %3.2f;...\n',p1,p2,p3);

   g1 = ( 1 - beta1 ) * trace( H11 * H11') * pS - e;
   g2 = ( 1 - beta2 ) * trace( H22 * H22') * pR - e;
   g3 = ( 1 - beta3 ) * trace( H33 * H33') * pO - e;
 
   lamdaS = lamda1 - tau1 * g1;
   lamdaR = lamda2 - tau2 * g2;
   lamdaO = lamda3 - tau3 * g3;

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

  lamda1 = lamdaS1;
  lamda2 = lamdaR1;
  lamda3 = lamdaO1;


    
   

   

     
end


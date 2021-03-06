function [ pS,pR,pO ] = WaterFilling_alg_3( H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,iter_num,P1_max,P2_max,P3_max,e,aL )
pS = zeros( 1, iter_num + 1 );
pR = zeros( 1, iter_num + 1 );
pO = zeros( 1, iter_num + 1 );

muS = zeros( 1, iter_num + 1 );
muR = zeros( 1, iter_num + 1 );
muO = zeros( 1, iter_num + 1 );

gammaS = zeros( 1, iter_num + 1 );
gammaR = zeros( 1, iter_num + 1 );
gammaO = zeros( 1, iter_num + 1 );

lamdaS = zeros( 1, iter_num + 1 );
lamdaR = zeros( 1, iter_num + 1 );
lamdaO = zeros( 1, iter_num + 1 );


%% Initialize 
p1 = 0.5;
p2 = 0.5;
p3 = 0.5;

pS( :,1 ) = p1;
pR( :,1 ) = p2;
pO( :,1 ) = p3;

lamda1 = aL;
lamda2 = aL;
lamda3 = aL;

lamdaS ( :,1 ) = lamda1;
lamdaR ( :,1 ) = lamda2;
lamdaO ( :,1 ) = lamda3;

tau1 = 0.0001;
tau2 = 0.0001;
tau3 = 0.0001;


%% obtain p1��p2 with lamda1,lamda2(gamma1,gamma2) fixed
% p1 = 1 \ ( omega1 + gamma1 ) - ( beta1 * norm( H12 * H12' ) * p2 + ( sigma_n )^2 ) \ ( beta1 * norm( H11 * H11' ) );
% p2 = 1 \ ( omega2 + gamma2 ) - ( beta2 * norm( H21 * H21' ) * p1 + ( sigma_n )^2 ) \ ( beta2 * norm( H22 * H22' ) );
% omega1 = ( beta1 * norm( H11 * H11' ) ) \ ( beta1 * norm( H11 * H11') * p1 + beta1 * norm( H12 * H12')  * p2 + ( sigma_n )^2 );
% omega2 = ( beta2 * norm( H22 * H22' ) ) \ ( beta2 * norm( H21 * H21') * p1 + beta2 * norm( H22 * H22')  * p2 + ( sigma_n )^2 );

%% obtain lamda1,lamda2 with p1��p2 fixed
% g1 = ( 1 - bata1 ) * norm( H11 * H11' ) * p1 - e1;
% g2 = ( 1 - bata2 ) * norm( H22 * H22' ) * p2 - e2;
% lamda1 = lamda1 - tau1 * g1;
% lamda2 = lamda2 - tau2 * g2;
% gamma1 = lamda1 * ( 1 - beta1 ) * norm( H11 * H11' );
% gamma2 = lamda2 * ( 1 - beta2 ) * norm( H22 * H22' );

%% The Interation Loop
for iter = 1 : iter_num
   
    fprintf( 'The %d th iteration ......',iter );
    fprintf( 'p1 = %3.2f; p2 = %3.2f;p3 = %3.2f; p4 = %3.2f;...',p1,p2,p3)
%% Update p1��p2 with lamda1,lamda2 fixed
  [ mu1,mu2,mu3 ] = get_mu3( p1,p2,p3,sigma_n,beta1,beta2,beta3,H11,H22,H33,H12,H13,H21,H23,H31,H32,P1_max,P2_max,P3_max );
  fprintf( 'mu1 = %5.4f; mu2 = %5.4f;mu3 = %5.4f; mu4 = %5.4f;...',mu1,mu2,mu3 );
  muS( :,iter ) = mu1;
  muR( :,iter ) = mu2;
  muO( :,iter ) = mu3;

  gammaS( :,iter ) = lamda1 * ( 1 - beta1 ) * trace( H11 * H11');
  gammaR( :,iter ) = lamda2 * ( 1 - beta2 ) * trace( H22 * H22');
  gammaO( :,iter ) = lamda3 * ( 1 - beta3 ) * trace( H33 * H33');

  A1 = ( muS( :,iter ) + gammaS( :,iter ) )^(-1) - ( beta1 * trace( H21 * H21') * pR( :,iter ) + beta1 * trace( H31 * H31') * pO( :,iter )  + sigma_n ) / ( beta1 * trace( H11 * H11') );
    if A1 <= 0
       pS( :,iter + 1 ) = 0;
    else 
       if A1 >= P1_max
          pS( :,iter + 1 ) = P1_max;
       else 
          pS( :,iter + 1 ) = A1;
       end
    end
    
    
 A2 = ( muR( :,iter ) + gammaR( :,iter ) )^(-1) - ( beta2 * trace( H12 * H12') * pS( :,iter + 1 ) + beta2 * trace( H32 * H32') * pO( :,iter )  + sigma_n ) / ( beta2 * trace( H22 * H22') );
    if A2 <= 0
       pR( :,iter + 1 ) = 0;
    else
       if A2 >= P2_max
          pR( :,iter + 1 ) = P2_max;
       else
          pR( :,iter + 1 ) = A2;
       end
    end
    
    
  A3 = ( muO( :,iter ) + gammaO( :,iter ) )^(-1) - ( beta3 * trace( H13 * H13') * pS( :,iter + 1 ) + beta3 * trace( H23 * H23') * pR( :,iter + 1 ) + + sigma_n ) / ( beta3 * trace( H33 * H33') );
    if A3 <= 0
       pO( :,iter + 1 ) = 0;
    else 
       if A3 >= P3_max
          pO( :,iter + 1 ) = P3_max;
       else 
          pO( :,iter + 1 ) = A3;
       end
    end
    p1 = pS( :,iter + 1 );
    p2 = pR( :,iter + 1 );
    p3 = pO( :,iter + 1 );
     %% Update lamda1,lamda2 with p1��p2 fixed
   fprintf( 'lamda1 = %5.4f; lamda2 = %5.4f;lamda3 = %5.4f; lamda4 = %5.4f;\n',lamda1,lamda2,lamda3);
   g1( :,iter ) = ( 1 - beta1 ) * trace( H11 * H11') * pS( :,iter ) - e;
   g2( :,iter ) = ( 1 - beta2 ) * trace( H22 * H22') * pR( :,iter ) - e;
   g3( :,iter ) = ( 1 - beta3 ) * trace( H33 * H33') * pO( :,iter ) - e;
   
   lamdaS( :,iter + 1 ) = lamdaS( :,iter ) - tau1 * g1( :,iter );
   lamdaR( :,iter + 1 ) = lamdaR( :,iter ) - tau2 * g2( :,iter );
   lamdaO( :,iter + 1 ) = lamdaO( :,iter ) - tau3 * g3( :,iter );
  
   B1 = lamdaS( :,iter ) - tau1 * g1( :,iter );     
       if B1 < 0
            lamdaS( :,iter + 1 ) = 0;
         else
            lamdaS( :,iter + 1 ) = B1;
       end
  
   B2 = lamdaR( :,iter ) - tau2 * g2( :,iter );
       if B2 < 0
            lamdaR( :,iter + 1 ) = 0;
         else
            lamdaR( :,iter + 1 ) = B2;
       end
    
   B3 = lamdaO( :,iter ) - tau3 * g3( :,iter );
       if B3 < 0
            lamdaO( :,iter + 1 ) = 0;
         else
            lamdaO( :,iter + 1 ) = B3;
       end
    

   lamda1 = lamdaS( :,iter + 1 );
   lamda2 = lamdaR( :,iter + 1 );
   lamda3 = lamdaO( :,iter + 1 );



end
end
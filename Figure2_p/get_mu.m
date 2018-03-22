function [ mu1,mu2,mu3,mu4 ] = get_mu( p1,p2,p3,p4,sigma_n,beta1,beta2,beta3,beta4,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,P1_max,P2_max,P3_max,P4_max )
%% Find lagrange multipliers, mu1 mu2 mu3 mu4
    syms mu1 mu2 mu3 mu4;
    
%% find the roots
equal1 = ( mu1 )^(-1) - ( beta1 * trace( H12 * H12') * p2 + beta1 * trace( H13 * H13') * p3 + beta1 * trace( H14 * H14') * p4 + sigma_n ) / ( beta1 * trace( H11 * H11') ) - P1_max;
roots1 = solve(equal1,'mu1');
roots1=double(roots1);
equal2 = ( mu2 )^(-1) - ( beta2 * trace( H21 * H21') * p1 + beta2 * trace( H23 * H23') * p3 + beta2 * trace( H24 * H24') * p4 + sigma_n ) / ( beta2 * trace( H22 * H22') ) - P2_max;
roots2 = solve(equal2,'mu2');
roots2=double(roots2);
equal3 = ( mu3 )^(-1) - ( beta3 * trace( H31 * H31') * p1 + beta3 * trace( H32 * H32') * p2 + beta3 * trace( H34 * H34') * p4 + sigma_n ) / ( beta3 * trace( H33 * H33') ) - P3_max;
roots3 = solve(equal3,'mu3');
roots3=double(roots3);
equal4 = ( mu4 )^(-1) - ( beta4 * trace( H41 * H41') * p1 + beta4 * trace( H42 * H42') * p2 + beta4 * trace( H43 * H43') * p3 + sigma_n ) / ( beta4 * trace( H44 * H44') ) - P4_max;
roots4 = solve(equal4,'mu4');
roots4=double(roots4);

%% get positive value mu according to KKT, Lagrange multiplier
for r = 1:length(roots1)
  if roots1(r) >0
     mu1 = roots1(r);
   end
end
for r = 1:length(roots2)
   if roots2(r) >0
      mu2 = roots2(r);
   end
end 
for r = 1:length(roots3)
  if roots3(r) >0
     mu3 = roots3(r);
   end
end
for r = 1:length(roots4)
   if roots4(r) >0
      mu4 = roots4(r);
   end
end 
end

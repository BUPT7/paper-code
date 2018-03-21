function [ mu1,mu2 ] = get_mu( p1,p2,sigma_n,beta1,beta2,H11,H22,H12,H21,P1_max,P2_max )
%% Find lagrange multipliers, mu1 mu2 mu3 mu4
    syms mu1 mu2;
    
%% find the roots
equal1 = ( mu1 )^(-1) - ( beta1 * trace( H21 * H21') * p2  +  sigma_n ) / ( beta1 * trace( H11 * H11') ) - P1_max;
roots1 = solve(equal1,'mu1');
roots1=double(roots1);
equal2 = ( mu2 )^(-1) - ( beta2 * trace( H12 * H12') * p1 + sigma_n ) / ( beta2 * trace( H22 * H22') ) - P2_max;
roots2 = solve(equal2,'mu2');
roots2=double(roots2);



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

end

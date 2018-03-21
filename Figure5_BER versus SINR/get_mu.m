function [ mu1,mu2 ] = get_mu( p1,p2,sigma_n,beta1,beta2,a,b,c,d,P1_max,P2_max )
%% Find lagrange multipliers, mu1 mu2 
    syms mu1 mu2 ;
    
%% find the roots
equal1 = ( mu1 )^(-1) - ( beta1 * a * p2 + sigma_n) / ( beta1 * c ) - P1_max;
roots1 = solve(equal1,'mu1');
roots1=double(roots1);
equal2 = ( mu2 )^(-1) - ( beta2 * b * p1 + sigma_n ) / ( beta2 * d ) - P2_max;
roots2 = solve(equal2,'mu2');
roots2=double(roots2);

%% get positive value mu according to KKT, Lagrange multiplier
mu1=0.001;
mu2=0.001;
for r = 1:length(roots1)
  if roots1(r) > 0
     mu1 = roots1(r);
   end
end
for r = 1:length(roots2)
   if roots2(r) > 0
      mu2 = roots2(r);
   end
end 
end

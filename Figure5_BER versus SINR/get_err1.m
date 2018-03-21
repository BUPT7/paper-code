function [ err_1,err_2,err_3,err ] = get_err1( p1_id,p2_id,p3_id,H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,c1,c2,c3,n_r )
%%  QPSK
x1_ = 2*c1-1;%  QPSK
x2_ = 2*c2-1;
x3_ = 2*c3-1;
% 
x1 = x1_(:,1) + 1i * x1_(:,2);
x2 = x2_(:,1) + 1i * x2_(:,2);
x3 = x3_(:,1) + 1i * x3_(:,2);

x1 = x1 / sqrt(2);
x2 = x2 / sqrt(2);
x3 = x3 / sqrt(2);

%%  ERR Calculation
y1 = sqrt(beta1*p1_id)*H11*x1 + sqrt(beta2*p2_id)*H21*x2 + sqrt(beta3*p3_id)*H31*x3 + n_r;             
y2 = sqrt(beta2*p2_id)*H22*x2 + sqrt(beta1*p1_id)*H12*x1 + sqrt(beta3*p3_id)*H32*x3 + n_r;             %Obtain the received signal y1 y2
y3 = sqrt(beta3*p3_id)*H33*x3 + sqrt(beta1*p1_id)*H13*x1 + sqrt(beta2*p2_id)*H23*x2 + n_r;

                                                              %Decode the received signal y1 y2

%% The Sum_ERR
err_1 = sqrt(2)*abs(y1-x1)>1;
err_2 = sqrt(2)*abs(y2-x2)>1;
err_3 = sqrt(2)*abs(y3-x3)>1;

err = sum(err_1) + sum(err_2)+ sum(err_3);
%err = sum(err_2);%Calculate the Total BER
end
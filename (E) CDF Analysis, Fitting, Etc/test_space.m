
AveWaitTime = 50; %min
DeadTime = 5;
NGammaSteps = 1;

x = 5.1:100;
DeadYVal = expcdf(DeadTime,AveWaitTime);

%test = (expcdf(x,AveWaitTime)-DeadYVal)./(1-DeadYVal);
test = (gamcdf(x,NGammaSteps,AveWaitTime)-DeadYVal)./(1-DeadYVal);
figure(47)
plot(x,test,'--')
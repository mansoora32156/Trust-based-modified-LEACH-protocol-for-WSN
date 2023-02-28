clc;
a1=0.9;
a2=10;
a3=4;
a4=0.5;
a5=10;
a6=4;
a7=0.9;
a8=15;
a9=4;
a10=0.9;
a11=10;
a12=2;
ACj=0:0.2:1;
len=length(ACj);
theta1 = zeros(1, len);
theta2 = zeros(1, len);
theta3 = zeros(1, len);
theta4 = zeros(1, len);
for i=1:len
theta1(i) = theta1(i)+(1-(a1/(1+exp(-a2*ACj(i)+a3))));
theta2(i) = theta2(i)+(1-(a4/(1+exp(-a5*ACj(i)+a6))));
theta3(i) = theta3(i)+(1-(a7/(1+exp(-a8*ACj(i)+a9))));
theta4(i) = theta4(i)+(1-(a10/(1+exp(-a11*ACj(i)+a12))));
disp(ACj);
disp(theta1);
disp(theta2);
disp(theta3);
disp(theta4);
plot(ACj,theta1,'b');
hold on;
plot(ACj,theta2,'g');
plot(ACj,theta3,'r');
plot(ACj,theta4,'k');
hold off;
title("Adaptive penalty coefficient under different parameters")
legend('a1=0.9;a2=10;a3=3','a1=0.5;a2=10;a3=3','a1=0.9;a2=15;a3=3','a1=0.9;a2=10;a3=2')
ylabel("Adaptive penalty coefficient");
xlabel("Proportion of abnormal behavior(theta)")
end


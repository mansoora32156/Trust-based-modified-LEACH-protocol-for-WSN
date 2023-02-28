clear all;
clc;
a1=0.9;
a2=10;
a3=4;
theta=zeros(1, 6);
    for ACj=0:0.2:1
        theta = 1-(a1/(1+exp(-a2*ACj+a3)));
        plot(theta,ACj); 
    end

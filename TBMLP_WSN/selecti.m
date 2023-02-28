clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Plot of Acj vs Theta                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

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
figure(1);
plot(ACj,theta1,'b');
hold on;
plot(ACj,theta2,'g');
plot(ACj,theta3,'r');
plot(ACj,theta4,'k');
hold off;
title("Adaptive penalty coefficient under different parameters")
legend('a1=0.9;a2=10;a3=3','a1=0.5;a2=10;a3=3','a1=0.9;a2=15;a3=3','a1=0.9;a2=10;a3=2')
ylabel("Adaptive penalty coefficient(theta)");
xlabel("Proportion of abnormal behavior(ACj)")
end


%function CT= comp_trust( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       Calculation of Direct trust & Indirect trust      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
m1=0.33;m2=0.33;m3=0.33; % n1,n2,n3 for CT
%n= 100;
%%%% Direct Trust
n = 100;
%r = 10;
pi = 3.14;
a1 = 2;
b1 = 1.25;
p_1 = [];
q_1 = [];
p_2 = [];
q_2 = [];
fr_1 = [];
fr_2 = [];
f_r = [];
err = [];
Eerr = [];
s1 = rand(1,100);
s1 = [s1];
s2 = rand(1,100);
s2 = [s2];
S = [];
r_c1 = [];
r_c2 = [];
r_c = [];
Mass_function = [];
m = [];
bel = [];
Belief_function = [];
i_t = [];
a_dt = [];
delta1 = rand(1,100);
delta1 = [delta1];
%disp("Delta1 = ")
%disp(delta1)
delta2 = rand(1,100);
delta2 = [delta2];
% disp("Delta2 = ")
% disp(delta2)
delta3 = rand(1,100);
delta3 = [delta3];
% disp("Delta3 = ")
% disp(delta3)
d_t1 = [];
d_t2 = [];
d_t = [];
for i = 1:n
    p_1(i) = i;
%     disp("Number of bits that are transmitted in the forward cycle=")
%     disp(p1(i))
    q_1(i) = n - p_1(i);
%     disp("Number of bits that are dropped in the forward cycle = ")
%     disp(q1(i))
    fr_1(i) = p_1(i)/(p_1(i)+q_1(i));
%     disp("The forwarding ratios in the forward cycle are = ")
%     disp(fr_1(i))         
end
for f1 = 1:length(fr_1)
    if fr_1(f1) < fr_1(f1+1:length(fr_1))
        delta1(f1+1:length(fr_1)) = delta1(f1) + (a1*(fr_1(f1) - fr_1(f1+1:length(fr_1))));
%         disp("The 'delta1' parameter values are = ")
%         disp(delta1(f1+1:length(fr1)));
    end
    %d_t1(f1+1:length(fr_1)) = fr_1(f1) * cos((pi/2) * delta1(f1+1:length(fr_1)));
    d_t1(f1+1:length(fr_1)) = fr_1(f1) * cos((pi/2) * delta1(f1+1:length(fr_1)));
    
    %d_t1(f1+1:length(fr_1)) = d_t1(f1) * cos((pi/2) * delta1(f1+1:length(fr_1)));
    
end

disp("The direct trust values in the forward cycle are = ");
disp(d_t1);    

for c = n:-1:1
    p_2(c) = c;
%     disp("Number of bits that are transmitted in the backward cycle = ")
%     disp(p2(j))
    q_2(c) = n - p_2(c);
%     disp("Number of bits that are dropped in the backward cycle = ")
%     disp(q2(j))
    fr_2(c) = p_2(c)/(p_2(c)+q_2(c));
%     disp("The forwarding ratios in the backward cycle are = ")
%     disp(fr_2(j))
end
for f2 = 1:length(fr_2)    
    if fr_2(f2) > fr_2(f2+1:length(fr_2))
        delta2(f2+1:length(fr_2)) = delta2(f2) + (b1*(fr_2(f2) - fr_2(f2+1:length(fr_2))));
%         disp("The 'delta2' parameter values are  = ")
%         disp(delta2(f2+1:length(fr2)));
    else
        delta3(f2+1:length(fr_2)) = delta3(f2);
%         disp("The 'delta3' parameter values are  = ")
%         disp(delta3(f2+1:length(fr_2)));  
%         disp(delta3(f2));
    end
    d_ta(f2+1:length(fr_2)) = fr_2(f2) * cos((pi/2) * delta2(f2+1:length(fr_2)));
    
    %d_ta(f2+1:length(fr_2)) = d_ta(f2) * cos((pi/2) * delta2(f2+1:length(fr_2))); 
    
%     disp("The direct trust values with respect to 'delta2' parameter = ");
%     disp(d_ta);

    d_tb(f2+1:length(fr_2)) = fr_2(f2) * cos((pi/2) * delta3(f2+1:length(fr_2)));

    %d_tb(f2+1:length(fr_2)) = d_tb(f2) * cos((pi/2) * delta3(f2+1:length(fr_2))); 
    
%     disp("The direct trust values with respect to 'delta3' parameter = ");
%     disp(d_tb);
end

%d_ta(f2+1:length(fr_2)) = fr_2(f2) * cos((pi/2) * delta2(f2+1:length(fr_2)));
disp("The direct trust values with respect to 'delta2' parameter = ");
disp(d_ta);
%d_tb(f2+1:length(fr_2)) = fr_2(f2) * cos((pi/2) * delta3(f2+1:length(fr_2)));
disp("The direct trust values with respect to 'delta3' parameter = ");
disp(d_tb);

d_t2 = [d_ta,d_tb];
disp("The direct trust values in the backward cycle are = ");
disp(d_t2);  

fr_2 = flip(fr_2,2);
%disp(fr2)
f_r = [fr_1,fr_2];
disp("The forwardng ratios are = ")
disp(f_r);

d_t = [d_t1,d_t2];
disp("The direct trust values are = ")
disp(d_t);

%%%% Recomendation Credibility

% Root Mean Square Error :

for e = 1:length(d_t)
    err(e) = abs(sqrt(mean((d_t(e)-d_t(e+1:length(d_t)).^2))));
    disp("The magnitude of the root mean square among the nodes are, Eerr = ");
    disp((err(e)));
end
% Eerr = [abs(err(x))];
% disp(Eerr)

% Similarity Parameter :

gamma = 0.35; %Tolerance Threshold
eta = 13;
psi = 11;


for l1 = 1:length((err))% 'l1' is the loop variable.
    if (err(l1)) < gamma
        for x = 1:length(s1)
            s1(x+1:length(s1)) = s1(x) + ((1-s1(x))/eta);
%             disp(s1(x+1:length(s1)))
        end
    end
end
disp("The similarity parameters when the root mean square error is less than the tolerance threshold is, s1 = ")
disp(s1);

for l2 = 1:length((err))
    if (err(l2)) >= gamma
        for y = 1:length(s2)
            s2(y+1:length(s2)) = s2(y) - (s2(y)/psi);
%             disp(s2(y+1:length(s2)))
        end
    end
end
disp("The similarity parameters when the root mean square error is less than the tolerance threshold is, s2 = ")
disp(s2);

s = [s1,s2];
disp("The overall similarity parameters are, s = ")
disp(s);

% Recomendation Credibility :

theta = 0.01;

for r1 = 1:length(s)
    if s(r1) > theta
        R1(r1) = 1 - (log(s(r1))/log(theta));
    end
end
r_c1 = [R1];
disp("The recomendation credibility values are, r_c1 = ")
disp(r_c1);

for r2 = 1:length(s)
    if s(r2) <= theta
        R2(r2) = 0;
    end
end
r_c2 = [R2];
disp("The recomendation credibility values are r_c2 = ")
disp(r_c2)
        
r_c = [r_c1,r_c2];
disp("The final set of all the recomendation credibility values are , r_c = ")
disp(r_c);

%%%% Indirect Trust

% Mass Function :

for r = 1:length(r_c)
    if r_c(r) > 0
        for d = 1:length(d_t)
            if d_t(d) >0
                for k = 1:n
                    Mass_function(k) = (r_c(randi(r)) * d_t(d))/(r_c(r) + r_c(randi(r)));
%                   m = [mass_function];
%                   disp("The mass function values for various nodes in the wireless sensor network is, m = ")
%                   disp(m);
                end
            end
        end
    end
end
m = [Mass_function];
disp("The mass function values for various nodes in the wireless sensor network is, Mass_function = ")
disp(m);

% Belief function :

for h = 1:length(m)
    bel(h) = m(h) * m(randi(h));
end
Belief_function = [bel];
disp("The values for the belief of one node over another, that is, the belief function values are, Belief_function = ");
disp(Belief_function);

% Indirect Trust :

for H = 1:length(m)
    for b = 1:length(Belief_function)
        i_t(H) = ((m(H) + m(randi(H))) - Belief_function(randi(b)));
    end
end
i_t;
disp("The indirect trust values are, i_t = ")
disp(i_t);

for a = 1:(length(d_t)-200)
    a_dt(a) = d_t(a);
end
a_dt;
disp("The final set of direct trust values are, a_dt = ")
disp(a_dt);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Energy trust             %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%E = EJ_trust(n);
%disp("E_J");
%disp(E);


    xm=300;  %Dimensions of x and y
    ym=300;
    zm=300;
    sink.x=0.5*xm;  %distance of base station from the network
    sink.y=0.5*ym;
    sink.z=0.5*zm;
    sink.x=100;
    sink.y=100;
    sink.z=100;
    %n=10;
    l=4000;
    Eelec= 50*10^(-9); %Eelec=50nJ
    Efs= 10*10^(-12);  %Efs=10pJ
    Emp=0.0013*10^(-12); %Emp=0.0013pJ
    E0=0.5; 
    d0=87; %d0=87m
    for i =1:n
        S(n+1).xd=sink.x;
        S(n+1).yd=sink.y;
        S(n+1).zd=sink.z;
        S(i).xd=rand(1,1)*xm;
        S(i).yd=rand(1,1)*ym;
        S(i).zd=rand(1,1)*zm;
        d(i)=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 + (S(i).zd-(S(n+1).zd) )^2 );
          %for j=i:n
            E_rcv(i)=  l*Eelec;
            if (d(i)< d0)
                E_s(i)= l*Eelec+l*Efs*(d(i)^2);
            elseif (d(i) >= d0)
                E_s(i)= l*Eelec+l*Emp*(d(i)^4);
   
            end 
          %end
     end
    disp("d:");
    disp(d);
        
    %d=0:20:100; %make changes for d accordingly (might use eucledian distance formula)
    
   
  
   
   
   
   ET=E_s;
   
   
   %disp("R_E");
   %disp(R_E);
   disp("Energy_Trust");
   disp(ET);
   %disp("Ej");
   %disp(Ej);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Comprehensive trust          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
ct = [];
C_T = [];
for N =1:n
    C_T(N)= m1*a_dt(N)+m2*i_t(N)+m3*ET(N);
end

C_T = [C_T];
disp("Comprehensive_Trust");
disp(C_T);
rmax=100;
r=1:rmax;
figure(4)
plot(r,C_T);
title("Comprehensive Trust(C_T) vs Round time(r)");
xlabel("Round Time");
ylabel("Comprehensive Trust");

figure(21);
nm_ct = [];
nmct = [];
n_mct = [];
nmc_t = [];
for n_m = 1:length(C_T)
    if C_T(n_m) > 0
        nm_ct(n_m) = C_T(n_m);
    end
end
nmct = [nm_ct];
n_mct = [nmct zeros(1, (100-length(nmct)))];
nmc_t = sort(n_mct);
plot(r,nmc_t);
title("Comprehensive Trust(C_T) vs Round time(r)for non-malicious nodes");
xlabel("Round time")
ylabel("Comprehensive Trust (non-malicious nodes)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           CT values against selective forwarding attack      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=100;
alpha=13.42;
beta =11.2;

for i =1:n1
    transfer_coeff(i)=alpha*rand(1,1)/beta*rand(1,1);
    
end
srth=mean(transfer_coeff);
srth


 CT1 = [];
 C_T1 = [];
 for i = 1:n1
     if transfer_coeff(i)<srth 
         
         C_T1(i)=C_T(i)- srth*0.35; %drop
         
     end
         if transfer_coeff(i)>srth
             
                 C_T1(i)=C_T(i)- srth*0.85; %drop
             
         end
 end
 
 
 CT1 = [C_T1];
 CT1=sort(C_T1,'descend');
 L=length(CT1);
 disp("CT values against selective forwarding attack");
 rmax=100;
 r=1:L;
 figure(21)
 plot(r,CT1);
 title("CT values against selective forwarding");
 xlabel("round time");
 ylabel("Comprehensive Trust");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             CT values against warmhole attack                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   

n1=100;
t1 = -4 + (4-(-4))*rand(1,1);
for i =1:n1
    t2_(i)=-4 + (4-(-4))*rand(1,1);
end

%RTT
for i =1:n1
    t3_(i)=t2_(i)-t1;
end

%Hopcount
disp("Hopcount");
h=1:1:n1;
for i =1:n1
    ts_(i)=t3_(i)/h(i);
end

%threshold for ts = avg(ts_(i))
tth=mean(ts_);

%Detection of warmhole link
%If the ts_(i)<tth (threshold) and also hc(i)=2 hopcount =2;
%then detect the route i as a warmhole link 
%else the route is safe
CT1 = [];
C_T1 = [];
for i =1:n1
    if ts_(i)<tth 
        C_T1(i)=C_T(i)-0.35;
    end
end
ts_ = [];
disp("ts_");
CT1 = [C_T1];
CT1=sort(C_T1,'descend');
L=length(CT1);
disp("CT values against warmhole attack");
%rmax=100;
r=1:L;
figure(19)
plot(r,CT1);
title("CT values against warmhole attack");
xlabel("round time");
ylabel("Comprehensive Trust");


mn_w = [];
mnwa = [];
m_n_wh_a = [];
for wct = 1:length(CT1)
    if CT1(wct) < 0
        mnwa(wct) = CT1(wct);
    end
end
mn_w = [mnwa];
m_n_wh_a = sort(mnwa,'descend');
disp("Comprehensive trust values of malicious nodes that are present due to Wormhole Attack are = ")
disp(mnwa);
figure(20);
plot(r,mn_w)
xlabel("Round Time")
ylabel("Comprehensive Trust values for the malicious nodes of the WSN under Wormhole attack")
title("Identifying the malicious nodes in the WSN under Wormhole attack");

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           LEACH                                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   



xm=300;       %diameters of sensor network
ym=300;
zm=300;
sink.x=0.5*xm;  %distance of base station from the network
sink.y=0.5*ym;
sink.z=0.5*zm;
sink.x=100;
sink.y=100;
sink.z=100;

n=100;
p=0.1;        %probibilty of a node to become cluster head
Eo=0.5;         %energy supplied to each node
ETX=50*0.000000001;          %transmiter energy per node
ERX=50*0.000000001;           %reciever energy per mode
Efs=10e-12;                %amplification energy when d is less than d0
Emp=0.0013e-12;                %amplification energy  when d is greater than d0
EDA=5*0.000000001;        %Data Aggregation Energy

rmax=1000;
do=sqrt(Efs/Emp);
Et=0;
Ct = rand(1,100);
Ct = [Ct];
 
for h=1:1
    S(n+1).xd=sink.x;
    S(n+1).yd=sink.y;
    S(n+1).zd=sink.z;
    Et=0;

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).zd=rand(1,1)*zm;
    ZR(i)=S(i).zd;
    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 + (S(i).zd-(S(n+1).zd) )^2 );
    S(i).distance=distance;
    S(i).G=0;
    S(i).type='N';
    S(i).E=Eo;
    Et=Et+S(i).E;
    S(i).C_T=C_T(i);
    Ct=Ct+S(i).C_T;
    figure(2)
    plot3(S(i).xd,S(i).yd,S(i).zd,'bd');
    text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
    title("Generation of Nodes");
    xlabel("Length of WSN");
    ylabel("Breadth of WSN");
    zlabel("Height of WSN");
    hold on  ; 
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(n+1).zd+1,num2str(n+1));
hold off ;

    Ct = rand(1,100);
    Ct = [Ct];
for c1 = 1:length(C_T)
    S(c1).xd=rand(1,1)*xm;
    XR(c1)=S(c1).xd;
    S(c1).yd=rand(1,1)*ym;
    YR(c1)=S(c1).yd;
    S(c1).zd=rand(1,1)*zm;
    ZR(c1)=S(c1).zd;
    distance=sqrt( (S(c1).xd-(S(n+1).xd) )^2 + (S(c1).yd-(S(n+1).yd) )^2 + (S(c1).zd-(S(n+1).zd) )^2 );
    S(c1).distance=distance;
    S(c1).G=0;
    S(c1).type='N';
    S(c1).E=Eo;
    Et=Et+S(c1).E;
    S(c1).C_T=C_T(c1);
    Ct=Ct+S(c1).C_T;    
    
        Ct(c)=Ct(c)+S(i).C_T(c);
    if  C_T(c1) > 0
        figure(11)
%                plot3(S(i).xd,S(i).yd,S(i).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
%                text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plot3(S(c1).xd,S(c1).yd,S(c1).zd,'bo');
        hold on;
        text(S(c1).xd+1,S(c1).yd-0.5,S(c1).zd+1,num2str(c1));
        title("Generation of Nodes that are not malicious");
        xlabel("Length of WSN");
        ylabel("Breadth of WSN");
        zlabel("Height of WSN");
        
    end   
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c1+1).zd+1,num2str(n+1));
hold off ;

for c2 = 1:length(C_T)
    S(c2).xd=rand(1,1)*xm;
    XR(c2)=S(c2).xd;
    S(c2).yd=rand(1,1)*ym;
    YR(c2)=S(c2).yd;
    S(c2).zd=rand(1,1)*zm;
    ZR(c2)=S(c2).zd;
    distance=sqrt( (S(c2).xd-(S(n+1).xd) )^2 + (S(c2).yd-(S(n+1).yd) )^2 + (S(c2).zd-(S(n+1).zd) )^2 );
    S(c2).distance=distance;
    S(c2).G=0;
    S(c2).type='N';
    S(c2).E=Eo;
    Et=Et+S(c2).E;
    S(c2).C_T=C_T(c2);
    Ct=Ct+S(c2).C_T;
    if C_T(c2) <= 0
        figure(12)
%             plot3(S(i).xd,S(i).yd,S(i).zd,'bo');
%             text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plot3(S(c2).xd,S(c2).yd,S(c2).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
        hold on ;
        text(S(c2).xd+1,S(c2).yd-0.5,S(c2).zd+1,num2str(c2));
        title("Generation of Nodes that are malicious");
        xlabel("Length of WSN");
        ylabel("Breadth of WSN");
        zlabel("Height of WSN");
        
    end
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c2+1).zd+1,num2str(n+1));
hold off ;

for c3 = 1:length(C_T)
    S(c3).xd=rand(1,1)*xm;
    XR(c3)=S(c3).xd;
    S(c3).yd=rand(1,1)*ym;
    YR(c3)=S(c3).yd;
    S(c3).zd=rand(1,1)*zm;
    ZR(c3)=S(c3).zd;
    distance=sqrt( (S(c3).xd-(S(n+1).xd) )^2 + (S(c3).yd-(S(n+1).yd) )^2 + (S(c3).zd-(S(n+1).zd) )^2 );
    S(c3).distance=distance;
    S(c3).G=0;
    S(c3).type='N';
    S(c3).E=Eo;
    Et=Et+S(c3).E;
    S(c3).C_T=C_T(c3);
    Ct=Ct+S(c3).C_T;
    figure(13)
    if C_T(c3) <= 0
        figure(13)
%             plot3(S(i).xd,S(i).yd,S(i).zd,'bo');
%             text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plot3(S(c3).xd,S(c3).yd,S(c3).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
        hold on ;
        text(S(c3).xd+1,S(c3).yd-0.5,S(c3).zd+1,num2str(c3));
    end
    if C_T(c3) > 0
        figure(13)
%                plot3(S(i).xd,S(i).yd,S(i).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
%                text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plot3(S(c3).xd,S(c3).yd,S(c3).zd,'bo');
        text(S(c3).xd+1,S(c3).yd-0.5,S(c3).zd+1,num2str(c3));        
    end
    title("Generation of malicious nodes and non-malicious nodes");
    xlabel("Length of WSN");
    ylabel("Breadth of WSN");
    zlabel("Height of WSN");
    
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c3+1).zd+1,num2str(n+1));
hold off ;

x_coordinates = [];
y_coordinates = [];
z_coordinates = [];
xyz_coordinates = [];
XYZ_COORDINATES =[];
for c4 = 1:length(C_T)
    S(c4).xd=rand(1,1)*xm;
    XR(c4)=S(c4).xd;
    S(c4).yd=rand(1,1)*ym;
    YR(c4)=S(c4).yd;
    S(c4).zd=rand(1,1)*zm;
    ZR(c4)=S(c4).zd;
    distance=sqrt( (S(c4).xd-(S(n+1).xd) )^2 + (S(c4).yd-(S(n+1).yd) )^2 + (S(c4).zd-(S(n+1).zd) )^2 );
    S(c4).distance=distance;
    S(c4).G=0;
    S(c4).type='N';
    S(c4).E=Eo;
    Et=Et+S(c4).E;
    S(c4).C_T=C_T(c4);
    Ct=Ct+S(c4).C_T;
    
    if C_T(c4) <= 0
%         figure(13)
%             plot3(S(i).xd,S(i).yd,S(i).zd,'bo');
%             text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        figure(14);
        plothandles = plot3(S(c4).xd,S(c4).yd,S(c4).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
        hold on;
        texthandles = text(S(c4).xd+1,S(c4).yd-0.5,S(c4).zd+1,num2str(c4));
       
    end
    if C_T(c4) > 0
%         figure(13)
%                plot3(S(i).xd,S(i).yd,S(i).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
%                text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plothandles = plot3(S(c4).xd,S(c4).yd,S(c4).zd,'bo');
        texthandles = text(S(c4).xd+1,S(c4).yd-0.5,S(c4).zd+1,num2str(c4));      
        x_coordinates(c4) = [S(c4).xd];
        y_coordinates(c4) = [S(c4).yd];
        z_coordinates(c4) = [S(c4).zd];
        hold off;
    end
    if C_T(c4) < 0 
        delete(plothandles);
        delete(texthandles);
    end
    title("Generation of Nodes by deleting the malicious nodes");
    xlabel("Length of WSN");
    ylabel("Breadth of WSN");
    zlabel("Height of WSN");
    hold on;
end

disp("X-Coordinates = ");
disp(x_coordinates);

disp("Y-Coordinates = ");
disp(y_coordinates);

disp("Z-Coordinates = ");
disp(z_coordinates);

xyz_coordinates = [x_coordinates,y_coordinates,z_coordinates];
disp("The X,Y & Z coordinates of the sensor nodes in the form of an one-dimensional array are = ")
disp(xyz_coordinates);
XYZ_COORDINATES = [nonzeros(x_coordinates),nonzeros(y_coordinates),nonzeros(z_coordinates)];
disp("The X,Y & Z coordinates of the sensor nodes are = ");
disp(XYZ_COORDINATES);

xc = [XYZ_COORDINATES(:,1)];
yc = [XYZ_COORDINATES(:,2)];
zc = [XYZ_COORDINATES(:,3)];

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c4+1).zd+1,num2str(n+1));
hold off ;


for c5 = 1:length(C_T)
    S(c5).xd=rand(1,1)*xm;
    XR(c5)=S(c5).xd;
    S(c5).yd=rand(1,1)*ym;
    YR(c5)=S(c5).yd;
    S(c5).zd=rand(1,1)*zm;
    ZR(c5)=S(c5).zd;
    distance=sqrt( (S(c5).xd-(S(n+1).xd) )^2 + (S(c5).yd-(S(n+1).yd) )^2 + (S(c5).zd-(S(n+1).zd) )^2 );
    S(c5).distance=distance;
    S(c5).G=0;
    S(c5).type='N';
    S(c5).E=Eo;
    Et=Et+S(c5).E;
    S(c5).C_T=C_T(c5);
    Ct=Ct+S(c5).C_T;
    
    
    if C_T(c5) < 0
%         figure(13)
%             plot3(S(i).xd,S(i).yd,S(i).zd,'bo');
%             text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        figure(15);
        plothandles1 = plot3(S(c5).xd,S(c5).yd,S(c5).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
        hold on;
        texthandles1 = text(S(c5).xd+1,S(c5).yd-0.5,S(c5).zd+1,num2str(c5));
        
    end
    if C_T(c5) <= 0 
         delete(plothandles1);
         delete(texthandles1);
    end
    if 0 < C_T(c5) < 1.25
%         figure(13)
%                plot3(S(i).xd,S(i).yd,S(i).zd,'ro','MarkerSize',8,'MarkerFaceColor','k');
%                text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
        plothandles1 = plot3(S(c5).xd,S(c5).yd,S(c5).zd,'bo');
        texthandles1 = text(S(c5).xd+1,S(c5).yd-0.5,S(c5).zd+1,num2str(c5));        
    end
    if C_T(c5) >= 1.25
        ct(c5) = C_T(c5);
        disp("The comprehensive trust values for the cluster heads are, ct = ");
        disp(ct(c5));        
        plothandles1 = plot3(S(c5).xd,S(c5).yd,S(c5).zd,'ro','MarkerSize',10,'MarkerFaceColor','y');
        texthandles1 = text(S(c5).xd+1,S(c5).yd-0.5,S(c5).zd+1,num2str(c5));
%         title("Generation of Nodes by deleting the malicious nodes and along with cluster heads");
        title("Generation of Nodes by deleting the malicious nodes and along with cluster heads");
        xlabel("Length of WSN");
        ylabel("Breadth of WSN");
        zlabel("Height of WSN");
        hold off;
    end
K = length(nonzeros(ct));
disp("The number of cluster heads are, K = ");
disp(K);
%     if CT(c5) < 0 
%         delete(plothandles1);
%         delete(texthandles1);
%     end
%     title("Generation of Nodes by deleting the malicious nodes and along with cluster heads");
%     xlabel("Length of WSN");
%     ylabel("Breadth of WSN");
%     zlabel("Height of WSN");
    hold on;
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c5+1).zd+1,num2str(n+1));
hold off ;

x_c = [];
coordinates = [];
for c_t = 1:1:length(nonzeros(nmc_t))
    if nmc_t(c_t) > 0
        S(c_t).xd=rand(1,1)*xm;
        XR(c_t)=S(c_t).xd;
        S(c_t).yd=rand(1,1)*ym;
        YR(c_t)=S(c_t).yd;
        S(c_t).zd=rand(1,1)*zm;
        ZR(c_t)=S(c_t).zd;
        x_c = [XR,YR,ZR];
%     disp(length(x_c));
%     disp(height(x_c));
%     disp(width(x_c));
        coordinates = reshape(x_c,[100,3]);
%     disp(height(coordinates));
%     disp(width(coordinates));
        distance=sqrt( (S(c_t).xd-(S(n+1).xd) )^2 + (S(c_t).yd-(S(n+1).yd) )^2 + (S(c_t).zd-(S(n+1).zd) )^2 );
        S(c_t).distance=distance;
        S(c_t).G=0;
        S(c_t).type='N';
        S(c_t).E=Eo;
        Et=Et+S(c_t).E;
        S(c_t).nmc_t=nmc_t(c_t);
        Ct=Ct+S(c_t).nmc_t;
            
%     if nmc_t(c_t) > 0
        idx = kmeans(XYZ_COORDINATES,K);
        figure(16);
        scatter3(XYZ_COORDINATES(:,1),XYZ_COORDINATES(:,2),XYZ_COORDINATES(:,3),100,idx,'filled');
        hold on;
    end
    title("Clustering of the WSN")
    xlabel("Length of WSN")
    ylabel("Breadth of WSN")
    zlabel("Height of WSN")
    
end

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(c_t+1).zd+1,num2str(n+1));
hold off ;


% c_h = [];
ch = [];
y_c = [];
cdt = [];
clh = [];
x_y_z = [XYZ_COORDINATES];
for ctx = 1:1:length(nonzeros(nmc_t))
    if nmc_t(ctx) > 0 
        S(ctx).xd=rand(1,1)*xm;
        XR(ctx)=S(ctx).xd;
        S(ctx).yd=rand(1,1)*ym;
        YR(ctx)=S(ctx).yd;
        S(ctx).zd=rand(1,1)*zm;
        ZR(ctx)=S(ctx).zd;
        y_c = [XR,YR,ZR];
%     disp(length(y_c));
%     disp(height(y_c));
%     disp(width(y_c));
        cdt = reshape(y_c,[100,3]);
%     disp(height(cdt));
%     disp(width(cdt));
        distance=sqrt( (S(ctx).xd-(S(n+1).xd) )^2 + (S(ctx).yd-(S(n+1).yd) )^2 + (S(ctx).zd-(S(n+1).zd) )^2 );
        S(ctx).distance=distance;
        S(ctx).G=0;
        S(ctx).type='N';
        S(ctx).E=Eo;
        Et=Et+S(ctx).E;
        S(ctx).nmc_t=nmc_t(ctx);
        Ct=Ct+S(ctx).nmc_t;
        
        
%     ch = [ct];
%     if 0 < C_T(ctx) < 1.25
        [idx,ch] = kmeans(x_y_z,K);
        figure(17);
        scatter3(x_y_z(:,1),x_y_z(:,2),x_y_z(:,3),50,idx,'filled');
%         plot3(ch(:,1),ch(:,2),ch(:,3),'p','MarkerSize',16,'MarkerFaceColor','k');
        hold on;
%         title("Clustering of the WSN with cluster head selection")
        
        title("Clustering of the WSN with cluster head selection");
        xlabel("Length of WSN");
        ylabel("Breadth of WSN");
        zlabel("Height of WSN");
        hold off;
    end
    
%     c_h = [ch zeros(1,(300-length(ch)))];
%     clh = reshape(c_h,[100,3]);
%     if C_T(ctx) >= 1.25

% plot3(ch(:,1),ch(:,2),ch(:,3),'p','MarkerSize',16,'MarkerFaceColor','k');

%     end
%     title("Clustering of the WSN with cluster head selection")
%     xlabel("Length of WSN")
%     ylabel("Breadth of WSN")
%     zlabel("Height of WSN")
    hold on;
end

plot3(ch(:,1),ch(:,2),ch(:,3),'p','MarkerSize',16,'MarkerFaceColor','k');

% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(ctx+1).zd+1,num2str(n+1));
hold off ;


ch = [];
c_h = [];
z_c = [];
c_dt = [];
c_dt_x = [];
c_dt_y = [];
c_dt_z = [];
c_d_t = [];
dist = [];
nm_coordinates = [x_y_z];
loc_nearest_ch_bs = [];
ctbd = [];
chlo = [];
P1 = [];
P2 = [];
PTS = [];
for cty = 1:1:length(nonzeros(nmc_t))
    if nmc_t(cty) > 0     
        S(cty).xd=rand(1,1)*xm;
        XR(cty)=S(cty).xd;
        S(cty).yd=rand(1,1)*ym;
        YR(cty)=S(cty).yd;
        S(cty).zd=rand(1,1)*zm;
        ZR(cty)=S(cty).zd;
        z_c = [XR,YR,ZR];
%     disp(length(y_c));
%     disp(height(y_c));
%     disp(width(y_c));
        c_dt = reshape(z_c,[100,3]);
%     disp(height(cdt));
%     disp(width(cdt));
        c_dt_x = c_dt(:,1);
        c_dt_y = c_dt(:,2);
        c_dt_z = c_dt(:,3);
        distance=sqrt( (S(cty).xd-(S(n+1).xd) )^2 + (S(cty).yd-(S(n+1).yd) )^2 + (S(cty).zd-(S(n+1).zd) )^2 );
        S(cty).distance=distance;
        S(cty).G=0;
        S(cty).type='N';
        S(cty).E=Eo;
        Et=Et+S(cty).E;
        S(cty).nmc_t=nmc_t(cty);
        Ct=Ct+S(cty).nmc_t;
        
%     ch = [ct];
        [idx,ch,sumd,D] = kmeans(nm_coordinates,K);
        figure(18);
        scatter3(nm_coordinates(:,1),nm_coordinates(:,2),nm_coordinates(:,3),50,idx,'filled');
        hold on;
    end
%     dist = reshape(D,[100,3]);
%     if C_T(cty) >= 1.25
%         plot3(S(cty).xd,S(cty).yd,S(cty).zd,'p','MarkerSize',16,'MarkerFaceColor','k');
%         hold on;    
% %         plot3(c_dt_x(:),c_dt_y(:),c_dt_z(:),'r-');
%     end
%     dist = reshape(D,[100,3]);
%     plot3(c_dt_x(:),c_dt_y(:),c_dt_z(:),'r-.');
    title("Clustering of the WSN along with paths")
    xlabel("Length of WSN");
    ylabel("Breadth of WSN");
    zlabel("Height of WSN");
    %hold on
end

plot3(ch(:,1),ch(:,2),ch(:,3),'p','MarkerSize',16,'MarkerFaceColor','k');
hold on;
line(ch(:,1),ch(:,2),ch(:,3),'color','red','LineWidth',1.5,'linestyle','--');


% BASE STATION :
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(cty+1).zd+1,num2str(n+1));
hold off;

bs = [S(n+1).xd,S(n+1).yd,S(n+1).zd];
for s0 = 1:size(ch,1)
%     disp(ch(s0,:));
    dist(s0) = [norm(bs-ch(s0,:))];
end
disp("The distances between the cluster-heads and the base station are = ")
disp(dist)
min_dist = min(dist);
disp("The distance between the nearest cluster-head to the base station is = ");
disp(min_dist);
    
for s1 = 1:size(ch,1)
    ctbd(s1) = [norm(bs-ch(s1,:))];
    if ctbd(s1) == min_dist
        chlo(s1,:) = [ch(s1,:)];
    end
end
disp("The location of the nearest cluster-head to the base station is =");
disp(chlo);
loc_nearest_ch_to_bs = [reshape(nonzeros(chlo),[1,3])];
% line(loc_nearest_ch_to_bs(:,1),loc_nearest_ch_to_bs(:,2),loc_nearest_ch_to_bs(:,3)'color','green','linestyle','-')

P1 = [loc_nearest_ch_to_bs(:,1),loc_nearest_ch_to_bs(:,2),loc_nearest_ch_to_bs(:,3)];
P2 = [S(n+1).xd,S(n+1).yd,S(n+1).zd];
PTS = [P1;P2];
line(PTS(:,1),PTS(:,2),PTS(:,3),'color','green','LineWidth',5,'linestyle','-');

hold on

data = cell(K,1);
locations = [];
for k = 1:K
    data{k} = [xc(idx==k),yc(idx==k),zc(idx==k)];
%     locations(k) = [xc(idx==k),yc(idx==k),zc(idx==k)];
end
locations = reshape(data,[1,k]);
disp("data = ")
disp(data);
disp("Locations = ")
disp(locations);

pts = double.empty(2,3,0);

for l0 = 1:length(locations)
%     for c_l = 1:size(ch,1)
    for l1 = 1:size(locations{l0})
        pts = [ch(l0,:);unique(locations{l0}(l1,:),'rows')];
        disp("pts = ")
        disp(pts);
        line(pts(:,1),pts(:,2),pts(:,3),'color','black','linestyle','-');
    end
end
% disp("pts = ")
% disp(pts);  

hold off ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countCHs=0;  %variable, counts the cluster head
cluster=1;  %cluster is initialized as 1
flag_first_dead=0; %flag tells the first node dead
flag_half_dead=0;  %flag tells the 10th node dead
flag_all_dead=0;  %flag tells all nodes dead
first_dead=0;
half_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
packets_TO_BS_per_round=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=0:1:rmax
    r;
    packets_TO_BS_per_round=0;
    %Operations for epochs
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end
    
    %hold off;
    
    %Number of dead nodes
    dead=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i=1:1:n
        %checking if there is a dead node
        if (S(i).C_T<=0.35)
            %plot(S(i).xd,S(i).yd,'red .');
            
            dead=dead+1;
            if (dead==1)
              if(flag_first_dead==0)
                 first_dead=r;
                 flag_first_dead=1;
              end
            end
            if(dead==0.5*n)
              if(flag_half_dead==0)
                  half_dead=r;
                  flag_half_dead=1;
              end
            end
            if(dead==n)
              if(flag_all_dead==0)
                  all_dead=r;
                  flag_all_dead=1;
              end
            end
            
            %hold on;
        end
        if S(i).C_T>0
            S(i).type='N';
        end
    end
    
        %plot(S(n+1).xd,S(n+1).yd,'x');
        STATISTICS.DEAD(h,r+1)=dead;
        STATISTICS.ALLIVE(h,r+1)=allive-dead;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   countCHs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).C_T>0.35)
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                %Election of Cluster Heads for normal nodes
                if ( temp_rand <= ( p/ ( 1 - p * mod(r,round(1/p)) )) )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    
                    S(i).type='C';
                    S(i).G=round(1/p)-1;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    C(cluster).zd=S(i).zd;
                    %plot(S(i).xd,S(i).yd,'k*');
                    
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 + (S(i).zd-(S(n+1).zd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    Z(cluster)=S(i).zd;
                    cluster=cluster+1;
                    
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                    end
                end
    
            end     
        end
    end
        
         STATISTICS.COUNTCHS(h,r+1)=countCHs;
    % or STATISTICS.COUNTCHS(h,r+1)=clster-1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Election of Associated Cluster Head for Normal Nodes
     for i=1:1:n
       if ( S(i).type=='N' && S(i).C_T>0.35)
        if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 + (S(i).zd-(S(n+1).zd) )^2 );
       min_dis_cluster=0;
         for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 + (S(i).zd-C(c).zd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %Calculating the culsterheads%
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
      
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH=packets_TO_CH+1;
       else 
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            packets_TO_BS_per_round=packets_TO_BS_per_round+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
       end
        S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   else
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 + (S(i).zd-S(n+1).zd)^2 );
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            packets_TO_BS_per_round=packets_TO_BS_per_round+1;
            
   end
  end
end
STATISTICS.PACKETS_TO_CH(h,r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(h,r+1)=packets_TO_BS;
STATISTICS.PACKETS_TO_BS_PER_ROUND(h,r+1)=packets_TO_BS_per_round;
STATISTICS.THROUGHPUT(h,r+1)=STATISTICS.PACKETS_TO_BS(h,r+1)+STATISTICS.PACKETS_TO_CH(h,r+1);

 En=0;
for i=1:n
    if S(i).E<=0
        continue;
    end
    En=En+S(i).E;
end
ENERGY(r+1)=En;
STATISTICS.ENERGY(h,r+1)=En;

end
first_dead_LEACH(h)=first_dead;
half_dead_LEACH(h)=half_dead;
all_dead_LEACH(h)=all_dead;

% cluster head display-------



% cluster head display-------
figure(3)
warning('OFF');
[Vx,Vy]=voronoi(X(:),Y(:));
plot(X,Y,'rd',Vx,Vy,'k-');
hold on;
title("Formation of cluster-heads");
xlabel("length of Wireless Sensor Networks(WSN)");
ylabel("height of Wireless Sensor Networks(WSN)");
voronoi(X,Y);
axis([10 xm 0 ym]);
hold off;

figure(5)
warning('OFF');
Vx=X(:);
Vy=Y(:);
Vz=Z(:);
plot3(X,Y,Z,'rd',Vx,Vy,Vz,'g-');
%plot(X,Y,'rd',vx,vy,'k--');
hold on;
title("Formation of path for transmission");
xlabel("Length of WSN");
ylabel("Breadth of WSN");
zlabel("Height of WSN");
%displaying the base station
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'d', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
hold off ;
%voronoi(X,Y,Z);
axis([10 xm 10 ym 10 zm]);
end

figure(6);
N=input("enter points ");
clf
hold on
plot(X,Y,'rd',Vx,Vy,'g-')
hold on
plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
[x,y]=ginput(N);
plot(x,y,'*')
title("Selection of Closest Node to BS ");
xlabel("Length of WSN");
ylabel("Breadth of WSN");
for i=1:1:N
    Dis(i)=sqrt((S(n+1).xd-x(i))^2+(S(n+1).yd-y(i))^2);
end
Dis=[Dis];
disp('Dis');
disp(Dis);
%Dis=sort(Dis)
min1=min(Dis);
i = find(Dis==min1);
%i=find(min1)
disp("closest Node to the BaseStation: ");
disp([x(i),y(i)]);
hold on
plot([x(i),S(n+1).xd],[y(i),S(n+1).yd],'k')
hold off


for r=0:rmax
    STATISTICS.DEAD(h+1,r+1)=sum(STATISTICS.DEAD(:,r+1))/h;
    STATISTICS.ALLIVE(h+1,r+1)=sum(STATISTICS.ALLIVE(:,r+1))/h;
    STATISTICS.PACKETS_TO_CH(h+1,r+1)=sum(STATISTICS.PACKETS_TO_CH(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS_PER_ROUND(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS_PER_ROUND(:,r+1))/h;
    STATISTICS.THROUGHPUT(h+1,r+1)=sum(STATISTICS.THROUGHPUT(:,r+1))/h;
    STATISTICS.COUNTCHS(h+1,r+1)=sum(STATISTICS.COUNTCHS(:,r+1))/h;
    STATISTICS.ENERGY(h+1,r+1)=sum(STATISTICS.ENERGY(:,r+1))/h;
end

first_dead=sum(first_dead_LEACH)/h;
half_dead=sum(half_dead_LEACH)/h;
all_dead=sum(all_dead_LEACH)/h;

%disp("S(i).CT");
%S.CT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0:rmax;
figure(7);
plot(r,STATISTICS.PACKETS_TO_BS(h+1,r+1));
title('pkts to BS')
xlabel("Round time(r)");
ylabel("pkts to BS")
figure(8);
plot(r,STATISTICS.PACKETS_TO_CH(h+1,r+1));
title('pkts to CH')
xlabel("Round time(r)");
ylabel("pkts to CH")
figure(9);
plot(r,STATISTICS.THROUGHPUT(h+1,r+1));
title('THROUGHPUT')
xlabel("Round time(r)");
ylabel("Throughput")
figure(10);
plot(r,STATISTICS.ENERGY(h+1,r+1));
title('Average Residual Energy') 
xlabel("Round time(r)");
ylabel("Residual Energy")
figure(22);
%Y1=STATISTICS.ENERGY/100;
%Y=1-Y1;
Y=100-STATISTICS.ENERGY;
plot(r,Y(h+1,r+1));
title('Average Energy Consumption') 
xlabel("Round time(r)");
ylabel("Energy consumed")   

  
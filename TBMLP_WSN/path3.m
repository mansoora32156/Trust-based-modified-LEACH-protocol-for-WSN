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
figure(13);
plot(ACj,theta1,'r');
hold on;
x=[-0.002801078,0.044817949,0.098039258,0.147058824,0.197478992,0.243697479,0.296918853,0.345938418,0.400560267,0.448179357,0.498599525,0.546218487,0.600840336,0.647058824,0.697478992,0.74789916,0.801120534,0.848739496,0.897759189,0.950980435,0.995798319];
y=[0.998333321,0.991,0.963,0.940,0.896,0.833,0.741,0.631,0.503,0.381,0.269,0.186,0.119,0.074,0.046,0.034,0.016,0.013,0.008,0.003,0.005];
plot(x,y,'b');
hold off;
legend('TBMLP', 'TBSEER');
title("Adaptive penalty coefficient under different parameters")
ylabel("Adaptive penalty coefficient(theta)");
xlabel("Proportion of abnormal behavior(ACj)")


%function CT= comp_trust( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Direct trust                 %%%%%%%%
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

for j = n:-1:1
    p_2(j) = j;
%     disp("Number of bits that are transmitted in the backward cycle = ")
%     disp(p2(j))
    q_2(j) = n - p_2(j);
%     disp("Number of bits that are dropped in the backward cycle = ")
%     disp(q2(j))
    fr_2(j) = p_2(j)/(p_2(j)+q_2(j));
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
            for k = 1:n
            Mass_function(k) = (r_c(randi(r)) * d_t(randi(d)))/(r_c(r) + r_c(randi(r)));
%             m = [mass_function];
%             disp("The mass function values for various nodes in the wireless sensor network is, m = ")
%             disp(m);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Indirect trust               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   


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
for j =1:n
    CT(j)= m1*a_dt(j)+m2*i_t(j)+m3*ET(j);
end


disp("Comprehensive_Trust");
disp(CT);
rmax=100;
r=1:rmax;
figure(4)
plot(r,CT);
title("Comprehensive Trust(CT) vs Round time(r)");
xlabel("rount time");
ylabel("Comprehensive Trust");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             CT values against warmhole attack               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
CT;
%CT1=zeros(1,n);
N=100;
t1 = -4 + (4-(-4))*rand(1,1);
for i =1:N
    t2_(i)=-4 + (4-(-4))*rand(1,1);
end

%RTT
for i =1:N
    t3_(i)=t2_(i)-t1;
end

%Hopcount
disp("Hopcount");
h=1:1:N;
for i =1:N
    ts_(i)=t3_(i)/h(i);
end

%threshold for ts = avg(ts_(i))
tth=mean(ts_);

%Detection of warmhole link
%If the ts_(i)<tth (threshold) and also hc(i)=2 hopcount =2;
%then detect the route i as a warmhole link 
%else the route is safe
for i =1:N
    if ts_(i)<tth 
        CT1(i)=CT(i)-0.35;
    end
end
disp("ts_");
CT1=sort(CT1,'descend');
L=length(CT1);
disp("CT values against warmhole attack");
%rmax=100;
r=1:L;
figure(11)
plot(r,CT1,'r');
title("CT values against warmhole attack");
xlabel("rount time");
ylabel("Comprehensive Trust");
hold on;
x=[-0.1439,11.952,24.192,35.999,47.807,59.903,71.855,83.808,96.0479,107.568];
y=[0.809, 0.569,0.503, 0.459, 0.427, 0.399, 0.383, 0.366, 0.353, 0.345];
plot(x,y,'b')
hold off;
legend('CT values for TBMLP','CT values for TBSEER');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           LEACH                             %%%%%%%%
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
Ct=0;
 
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
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).E=Eo;
    Et=Et+S(i).E;
    S(i).CT=CT(i);
    Ct=Ct+S(i).CT;
    figure(2)
      plot3(S(i).xd,S(i).yd,S(i).zd,'bo');
      text(S(i).xd+1,S(i).yd-0.5,S(i).zd+1,num2str(i));
      title("Generation of Nodes");
      xlabel("Length of WSN");
      ylabel("Breadth of WSN");
      zlabel("Height of WSN");
      hold on;
   
end

plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd+1,S(i).zd+1,num2str(n+1));
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
        if (S(i).CT<=0.35)
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
        if S(i).CT>0.35
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
        if(S(i).CT>0.35)
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
       if ( S(i).type=='N' && S(i).CT>0.35 )
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
%[vx,vy]=voronoi(X(:),Y(:));
figure(5)
warning('OFF');
vx=X(:);
vy=Y(:);
vz=Z(:);
plot3(X,Y,Z,'rd',vx,vy,vz,'g-');
%plot(X,Y,'rd',vx,vy,'k--');
hold on;
title("Formation of path for transmission");
xlabel("Length of WSN");
ylabel("Breadth of WSN");
zlabel("Height of WSN");
%displaying the base station
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
hold off ;
%voronoi(X,Y,Z);
axis([10 xm 10 ym 10 zm]);

end 
figure(6);
N=input("enter points ")
clf
hold on
plot(X,Y,'rd',vx,vy,'g-')
hold on
plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
[x,y]=ginput(N)
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
V=STATISTICS.ENERGY;
B=V(1,1);
B1=V/B;
plot(r,B1(h+1,r+1));
title('Average Residual Energy') 
xlabel("Round time(r)");
ylabel("Residual Energy");
figure(15);
V1=100-STATISTICS.ENERGY;
B1=V1(1,1);
B2=V1/B1;
plot(r,B2(h+1,r+1));
title('Average Energy Consumption') 
xlabel("Round time(r)");
ylabel("Energy consumed");
figure(12);
%Y1=STATISTICS.ENERGY/100;
%Y=1-Y1;
r1=0:500;
Y=100-STATISTICS.ENERGY;
Ymax=Y(1,1001);
Y1=Y/Ymax;
plot(r1,Y1(h+1,r1+1),'r');
title('Average Energy Consumption') 
xlabel("Round time(r)");
ylabel("Energy consumed");
hold on;
x=[18.232,59.607,79.943,100.981,121.318,140.252,161.290,181.626,201.963,218.092,238.429,258.064,281.907,298.036,319.074,338.008,362.552,379.382,399.018,417.251,438.990,459.326,476.858,500];
y=[0.238,0.516,0.715,0.794,0.874,0.993,1.152,1.391,1.550,1.748,2.066,2.066,2.225,2.384,2.543,2.702,2.821,3.060,3.179,3.219,3.418,3.656,3.855,3.974];
plot(x,y,'b');
legend('TBMLP', 'TBSEER')
hold off;
    
    
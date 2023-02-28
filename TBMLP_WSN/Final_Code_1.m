clear all;
clc;
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
%%%%                 Calculation of Direct trust                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    % 1) Calculation of Direct Trust Value
    % parameters
    % gamma (g) =0.5
    % lambda (l) =0.5
    m1=0.60;m2=0.30;m3=0.10; % n1,n2,n3 for CT
    n= 100;
    g=0.5;
    l=0.5;
    a1=0.9;
    a2=10;
    a3=4;
    ACj=0:0.2:n/5;
    len=length(ACj);
    theta = zeros(1, len);
    for i=1:len
        theta(i) = theta(i)+(1-(a1/(1+exp(-a2*ACj(i)+a3))));
    end

    for i= 1:n
        for j =i:n
            messagej=4400;
            receive_messagej=400;
            send_messagej=4000;
            rejectionj=send_messagej-receive_messagej;
            un_sendj=messagej-send_messagej;
            
            Rj=(theta*receive_messagej-rejectionj)/messagej;
            Sj=(theta*send_messagej-un_sendj)/messagej;
        end
    end
    
    DT=zeros(1,n).*0.5;
    %DT=abs(DT);
    HT=ones(1,n);
    IT=zeros(1,n);
    
    for i=1:n
            HT(i+1)=(l*(HT(i)+DT(i)));
            DT(i) = g*HT(i)+(1-g)*(Rj(i)+Sj(i));
    end
    disp("HT");
    disp(HT);
    DT=abs(DT);
    disp("Direct_Trust")
    disp(DT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Indirect trust               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   

    Bh = randsample(n,n)';
    for u = 1 : length(Bh)
    for q = u: n
        IT(u) = IT(u)+(1/q)*(DT(i)*DT(j));
    end
    end
    disp('Indirect_Trust');
    disp(IT);
    
    
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
    CT(j)= m1*DT(j)+m2*IT(j)+m3*ET(j);
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
%[vx,vy]=voronoi(X(:),Y(:));
figure(5)
warning('OFF');
vx=X(:);
vy=Y(:);
vz=Z(:);

plot3(X,Y,Z,'rd',vx,vy,vz,'g-');
%plot(X,Y,'rd',vx,vy,'k--');
hold on;
title("Formation of cluster-heads");
xlabel("Length of WSN");
ylabel("Breadth of WSN");
zlabel("Height of WSN");
%displaying the base station
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
hold off ;
%voronoi(X,Y,Z);
%voronoi(X,Y,Z);
axis([10 xm 10 ym 10 zm]);

end 
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
figure(6);
plot(r,STATISTICS.PACKETS_TO_BS(h+1,r+1));
title('pkts to BS')
xlabel("Round time(r)");
ylabel("pkts to BS")
figure(7);
plot(r,STATISTICS.PACKETS_TO_CH(h+1,r+1));
title('pkts to CH')
xlabel("Round time(r)");
ylabel("pkts to CH")
figure(8);
plot(r,STATISTICS.THROUGHPUT(h+1,r+1));
title('THROUGHPUT')
xlabel("Round time(r)");
ylabel("Throughput")
figure(9);
plot(r,STATISTICS.ENERGY(h+1,r+1));
title('Average Residual Energy') 
xlabel("Round time(r)");
ylabel("Residual Energy")
    
    
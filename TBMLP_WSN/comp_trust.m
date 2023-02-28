clear all;
clc;
%function CT= comp_trust( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Direct trust                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    % 1) Calculation of Direct Trust Value
    % parameters
    % gamma (g) =0.5
    % lambda (l) =0.5
    m1=0.66;m2=0.33;m3=0.00; % n1,n2,n3 for CT
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
    
    DT=ones(1,n).*0.5;
    HT=ones(1,n);
    IT=zeros(1,n);
    
    for i=1:n
            HT(i+1)=(l*(HT(i)+DT(i)));
            DT(i) = g*HT(i)+(1-g)*(Rj(i)+Sj(i));
    end
    disp("HT");
    disp(HT); 
    disp("DT");
    DT=abs(DT);
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
    disp('IT');
    disp(IT);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Energy trust             %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%E = EJ_trust(n);
%disp("E_J");
%disp(E);


    xm=300;  %Dimensions of x and y
    ym=300;
    sink.x=0.5*xm;  %distance of base station from the network
    sink.y=0.5*ym;
    sink.x=100;
    sink.y=75;
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
        S(i).xd=rand(1,1)*xm;
        S(i).yd=rand(1,1)*ym;
        d(i)=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
    end
    disp("d:");
    disp(d);
        
    %d=0:20:100; %make changes for d accordingly (might use eucledian distance formula)
    
   for j=i:n
   E_rcv(j)=  l*Eelec;
   if (d(j)< d0)
      E_s(j)= l*Eelec+l*Efs*(d(j)^2);
   elseif (d(j) >= d0)
      E_s(j)= l*Eelec+l*Emp*(d(j)^4);
   
   end 
   
   R_E(j)=E0-E_rcv(j)-E_s(j);
   
   E(j)=R_E(j)/E0;
   end
   disp("E_rcv");
   disp(E_rcv);
   disp("E_s");
   disp(E_s);
   disp("R_E");
   disp(R_E);
   disp("E_j");
   disp(E);
   %disp("Ej");
   %disp(Ej);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Calculation of Comprehensive trust          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for j =1:n
    CT(j)= m1*DT(j)+m2*IT(j)+m3*E(j);
end
disp("CT");
disp(CT);
    
            


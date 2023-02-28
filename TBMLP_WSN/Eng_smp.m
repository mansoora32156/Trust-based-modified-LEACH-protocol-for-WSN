clc;
% 3) Calculation of Energy Trust Value
    % parameters
    % radio frequency energy consuption of nodes (Eelec) = 50 nJ/ bit
    % length of message (l) = 4000 bits
    xm=300;  %Dimensions of x and y
    ym=300;
    sink.x=0.5*xm;  %distance of base station from the network
    sink.y=0.5*ym;
    sink.x=100;
    sink.y=75;
    n=6;
    l=4000;
    Eelec= 50*10^(-9); %Eelec=50nJ
    Efs= 10*10^(-12);  %Efs=10pJ
    Emp=0.0013*10^(-12); %Emp=0.0013pJ
    E0=0.5; 
    d0=87; %d0=87m
    
        
    d=0:20:100; %make changes for d accordingly (might use eucledian distance formula)
    
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
   
  %#mansoora pls plot tis xD
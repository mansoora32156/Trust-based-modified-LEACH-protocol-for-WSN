clc;
clear all;
close all;
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

gamma = 0.25; %Tolerance Threshold
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

% Indirect Trust :

for H = 1:length(m)
    for b = 1:length(Belief_function)
        i_t(H) = ((m(H) + m(randi(H))) - Belief_function(randi(b)));
    end
end
i_t = [i_t];
disp("The indirect trust values are, i_t = ")
disp(i_t);

for a = 1:(length(d_t)-200)
    a_dt(a) = d_t(a);
end
a_dt = [a_dt];
disp("The final set of direct trust values are, a_dt = ")
disp(a_dt);
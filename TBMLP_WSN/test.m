t = 100;
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
for i = 1:t
    p_1(i) = i;
%     disp("Number of bits that are transmitted in the forward cycle=")
%     disp(p1(i))
    q_1(i) = t - p_1(i);
%     disp("Number of bits that are dropped in the forward cycle = ")
%     disp(q1(i))
    fr_1(i) = p_1(i)/(p_1(i)+q_1(i));
    disp("The forwarding ratios in the forward cycle are = ")
    disp(fr_1(i))         
end
for f1 = 1:length(fr_1)
    if fr_1(f1) < fr_1(f1+1:length(fr_1))
        delta1(f1+1:length(fr_1)) = delta1(f1) + (a1*(fr_1(f1) - fr_1(f1+1:length(fr_1))));
%         disp("The 'delta1' parameter values are = ")
%         disp(delta1(f1+1:length(fr1)));
    end
    %d_t1(f1+1:length(fr_1)) = fr_1(f1) * cos((pi/2) * delta1(f1+1:length(fr_1)));
    d_t1(f1+1:length(fr_1)) = fr_1(f1) * cos((pi/2) * delta1(f1+1:length(fr_1)));
end

disp("The direct trust values in the forward cycle are = ");
disp(d_t1);    

for j = t:-1:1
    p_2(j) = j;
%     disp("Number of bits that are transmitted in the backward cycle = ")
%     disp(p2(j))
    q_2(j) = t - p_2(j);
%     disp("Number of bits that are dropped in the backward cycle = ")
%     disp(q2(j))
    fr_2(j) = p_2(j)/(p_2(j)+q_2(j));
    disp("The forwarding ratios in the backward cycle are = ")
    disp(fr_2(j))
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
%     disp("The direct trust values with respect to 'delta2' parameter = ");
%     disp(d_ta);
    d_tb(f2+1:length(fr_2)) = fr_2(f2) * cos((pi/2) * delta3(f2+1:length(fr_2)));
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

% % for f = 1:length(fr)
% %     if fr(f) < fr(f+1:length(fr))
% %         delta1(f+1:length(fr)) = delta1(f:length(fr)) + (a*(fr(f:length(fr)) - fr(f+1:length(fr))));
% %         disp("The 'delta1' parameter values are = ")
% %         disp(delta1(f+1:length(fr)));
% %     end 
% %     if fr(f) > fr(f+1:length(fr))
% %         delta2(f+1:length(fr)) = delta2(f) + (b*(fr(f) - fr(f+1:length(fr))));
% %         disp("The 'delta2' parameter values are  = ")
% %         disp(delta2(f+1:length(fr)));
% %     else
% %         delta3(f+1:length(fr)) = delta3(f);
% %         disp("The 'delta3' parameter values are  = ")
% %         disp(delta3(f+1:length(fr)));  
% %     end
% %     d_t(f+1:length(fr)) = fr(f) * cos((pi/2) * delta(f+1:length(fr)));
% % end
% d_t = [d_t1,d_t2];
% disp("The direct trust values are = ")
% disp(d_t);
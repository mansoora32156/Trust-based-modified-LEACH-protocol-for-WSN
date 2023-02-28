N=10;
t1 = -4 + (4-(-4))*rand(1,1);
for i =1:1:N
    t2_(i)= -4 + (4-(-4))*rand(1,1)
end
for i =1:N
    t3_(i)=t2_(i)-t1
end
h=1:1:N;
for i =1:N
    ts_(i)=t3_(i)/h(i)
end

tth=mean(ts_)
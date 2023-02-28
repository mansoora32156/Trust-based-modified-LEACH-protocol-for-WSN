a=linspace(0,2*pi);
b=sin(a);
n=input("enter points ")
clf
hold on
plot(a,b)
[x,y]=ginput(n)
plot(x,y,'*')

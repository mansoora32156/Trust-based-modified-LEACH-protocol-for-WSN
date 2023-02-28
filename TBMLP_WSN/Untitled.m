clc;
clear all;
close all; 
warning off;
disp('----------------------------------'); 


disp('LEACH-SAGA PROTOCOL FOR ENERGY BALANCE OF WIRELESS SENSOR NETWORK IN CLUSTERING ROUTING');

% pause (5); 
% detail info; 
% pause (5);
pop=input('Enter the population Value:--'); 
gen=input('Enter the Generation Value:--');
nodes=input('Enter the Number of Nodes:--'); 
net=input('Enter the Size Net:--'); 
packet=input('Enter the Packet Size:--');


no_of_nodes=nodes;%300; 
size_net=net; %300; 
packet_size=packet;%60; 
energy_con=6; 
xloc=randsrc(1, no_of_nodes, 1:size_net); 
yloc=randsrc(1, no_of_nodes, 1:size_net);
im=yloc;
kmin=9; % minimum number of classes 
kmax=9; % maximum number of classes
no_of_pop=pop;
final_database=-1*ones (no_of_pop, (kmax+ 2) +2);
for k3=1:no_of_pop

clst_size=9; 
centr_val=randsrc(1, clst_size,1:size_net); 
centr_vall=randsrc(1, clst_size, 1:size_net); 
[fitnessdb outclsdata lendsx]=FITNESS_PROCESS(xloc,yloc, centr_val, centr_vall);

final_database(k3,[1: (length(centr_val)*2) end-1:end]) =[centr_val centr_vall clst_size fitnessdb]; 

end
data_gen_process=final_database;

pc=0.8; 
pm=0.02;

no_of_gen=gen;

final_database2=-1*ones (no_of_pop, (kmax+2)+2); 
res_fin=ones (1, no_of_gen) *no_of_nodes; 
res_fin2=(1:(no_of_nodes*packet_size)) *packet_size*5; 
res_fin3=(1: (no_of_nodes*energy_con));
for iter=1:no_of_gen

    final_gaout=cross_over_process (data_gen_process, no_of_pop, pc, pm); % cross over & mutation process process going on

    for k4=1:no_of_pop

    locnew=9;
    centr_valnew=final_gaout (k4,1:locnew); 
    centr_valnewl=final_gaout (k4,locnew+1:(kmax*2));
    [fitnessdb outclsdata lendsx]=FITNESS_PROCESS(xloc, yloc, centr_valnew, centr_valnewl); 
    final_database2 (k4, [1: (length(centr_valnew)*2) end-1:end] )=[centr_valnew centr_valnewl length(centr_valnew) fitnessdb];
    datalen(k4)=lendsx;

end

data_gen_process=final_database2; 
final_database2=-1*ones (no_of_pop, (kmax+2)+2);
datax_deb(iter)=mean (datalen);

end 
datax_deb=sort(datax_deb); 
res_fin(end-(no_of_gen/2)-9:end-(no_of_gen/2))=datax_deb (1:10); 
res_fin(1:end-(no_of_gen/2) -9)=0; 
centr_valnewloc=data_gen_process(1, end-1); 
centr_valnewfinal=data_gen_process(1,1:centr_valnewloc); 
centr_valnewfinall=data_gen_process(1,centr_valnewloc+1:(kmax*2)); 
loc=find (res_fin==no_of_nodes);

res_fin2(1:loc (1))=res_fin2 (1:loc(1)); 
res_fin2(loc(1):loc(end))=no_of_nodes*packet_size; 
res_fin3(1:loc(1))=res_fin3 (1:loc(1)); 
res_fin3(loc(1):loc(end)) = (no_of_nodes/10) *energy_con;

for k4=1:length(centr_valnewfinal)
    fdx=centr_valnewfinal(k4); 
    fdy=centr_valnewfinal1(k4); 

    for k5=1:length(xloc)

        distfin(k4, k5)=sqrt((xloc (k5) -fdx).^2+ (yloc (k5)-fdy).^2); 
    end

end

[minval minloc]=min(distfin); 
outcls=minloc; figure, plot(xloc, yloc, 'bo', 'linewidth', 2); 
hold on; 
for k3=1:no_of_nodes
    text(xloc(k3) +0.3, yloc(k3), num2str(outcls(k3)), 'color','m','fontsize',20) 
end 
hold on, plot(centr_valnewfinal, centr_valnewfinal1, 'co', 'linewidth',3); 
hold on, plot(size_net/2, size_net/2, 'bo', 'linewidth',5); 
for k3=1:9
    text(centr_valnewfinal(k3) +0.3, centr_valnewfinall(k3), num2str((k3)), 'color', 'y','fontsize', 20) 
    hold on, plot([size_net/2 centr_valnewfinal(k3)], [size_net/2 centr_valnewfinal1(k3)], 'r:', 'linewidth', 2);

end
for k3=1:max(outcls)
    locmd=find (outcls==k3);

        clrrng=[rand rand rand] ;


    for k4=1:length(locmd)

        hold on, plot([xloc(locmd(k4)) centr_valnewfinal(k3)], [yloc(locmd(k4)) centr_valnewfinall(k3)], 'color', clrrng, 'linewidth', 2); 
    end
end
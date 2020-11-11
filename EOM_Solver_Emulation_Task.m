% Emulation Task
clear
clc
close all


tic


rng (17,'philox')

num_mass=10;
x_init = unifrnd(0,10,[num_mass 1]);
y_init = unifrnd(0,10,[num_mass 1]);
DT = delaunayTriangulation(x_init,y_init);

dt=delaunay(x_init,y_init);
figure
triplot(dt,x_init,y_init);
hold on
springs=[];
mass_locations=DT.Points; % initial location of the masses: the row corresponds to the mass index
springs(:,1:2)=DT.edges; % firs column is the first mass index and the second one is the second mass index for all the connected masses
xlim([0,10]);
ylim([0,10]);


for spring_counter=1:length(springs)
   
    % Third column of [springs] matrix will be the free length of each spirng connecting two corresponing masses 
   
   springs(spring_counter,3)=pdist([mass_locations(springs(spring_counter,1),1),mass_locations(springs(spring_counter,1),2);...
       mass_locations(springs(spring_counter,2),1),mass_locations(springs(spring_counter,2),2)],'euclidean'); 
    
end


range_rand=log10([100 200]);


springs(:,4)=unifrnd(0,100,[spring_counter 1]); %k1 values of each connection
springs(:,5)=10.^(range_rand(1) + (range_rand(2)-range_rand(1)).*rand(spring_counter,1)); %k3 values of each connection
springs(:,6)=unifrnd(0,100,[spring_counter 1]); %d1 values of each connection
springs(:,7)=10.^(range_rand(1) + (range_rand(2)-range_rand(1)).*rand(spring_counter,1)); %d3 values of each connection




all_springs=[springs;[springs(:,2),springs(:,1),springs(:,3),springs(:,4),springs(:,5),springs(:,6),springs(:,7)]];
all_springs_sorted=sortrows(all_springs,1);




%%


left_mass=find(mass_locations(:,1)==min(mass_locations(:,1)));
right_mass=find(mass_locations(:,1)==max(mass_locations(:,1)));
for mass_counter=1:num_mass
    if mass_counter==left_mass
        plot(mass_locations(mass_counter,1), mass_locations(mass_counter,2), '.r', 'MarkerSize',20)
    elseif mass_counter==right_mass
        plot(mass_locations(mass_counter,1), mass_locations(mass_counter,2), '.r', 'MarkerSize',20)
    else
        plot(mass_locations(mass_counter,1), mass_locations(mass_counter,2), '.b', 'MarkerSize',20)    
    end
end


initial_conditions=zeros(1,num_mass*4);
initial_conditions(1:4:end)=mass_locations(:,1);
initial_conditions(3:4:end)=mass_locations(:,2);

time_span=0:.001:100;

num_inputs=(.4)*num_mass;
%num_inputs=2;

mass_vec=1:num_mass;
active_masses=mass_vec(mass_vec~=left_mass);
active_masses=active_masses(active_masses~=right_mass);
input_nodes=randi([1,num_mass-2],num_inputs,1);
input_nodes=active_masses(input_nodes);
input_nodes_modified=zeros(num_mass,2);
input_nodes_modified(:,1)=1:num_mass;
input_nodes_modified(input_nodes,2)=(rand([num_inputs,1])*2)-1;


k1_vec=all_springs_sorted(:,4);
k3_vec=all_springs_sorted(:,5);
d1_vec=all_springs_sorted(:,6);
d3_vec=all_springs_sorted(:,7);


toc


tic

amp_target=2.5;
[T,P] = ode45(@(t,p) EOM_RECURRENT_SYSTEM_NODAMPING_VECTORIZED_FAST(t,p,all_springs_sorted,input_nodes_modified,num_mass,left_mass,right_mass,amp_target), time_span, initial_conditions);


toc

tic

figure
plot(T,P(:,1:4:end))
title('X-dir')

figure
plot(T,P(:,3:4:end))
title('Y-dir')
%%


%calculating x of each spring

mass_1=all_springs(1:size(all_springs,1)/2,1);
mass_2=all_springs(1:size(all_springs,1)/2,2);

spring_length=sqrt((P(1:length(T),4*mass_1(:)-3)-P(1:length(T),4*mass_2(:)-3)).^2+(P(1:length(T),4*mass_1(:)-1)-P(1:length(T),4*mass_2(:)-1)).^2);

figure
plot(T,normalize(spring_length(:,:)))   %plotting the normalized version around 0 and std=1.


toc



%%
%defining the state matrix for the learning phase
L=spring_length;
L=[ones(length(T(50000:75000)),1),L(50000:75000,:)];

%%

% dynamic system with memory (1)
% f1=2.11;
% f2=3.73;
% f3=4.33;
% f4=5.98;
% 
% 
% y_0=0;
% y=[y_0];
% time_span=0:.001:100;
% u_input=sin(2*pi*f1*time_span).*sin(2*pi*f2*time_span).*sin(2*pi*f3*time_span);
% for i=1:length(time_span)
%     
%    y_t=y(end)+2*u_input(i); 
%    y=[y;y_t]; 
%        
% end
% 
% target_function=(.1)*y(2:end);
% 
% Target=target_function(50000:75000);
% 
% 
% w_out=regress(Target,L);
% 
% for time_counter=1:length(T)
%     
%     output_signal(time_counter)=spring_length(time_counter,:)*w_out(2:end)+w_out(1);
%     
% end

%%
% mass-spring-damper system with memory (2)

% time_vec = 0:0.001:100;
% IC = [0 0];
% [t, y] = ode45(@msd_system_ode, time_vec, IC);
% target_function=y(:,1);
% 
% Target=target_function(50000:75000);
% 
% 
% w_out=regress(Target,L);
% 
% for time_counter=1:length(T)
%     
%     output_signal(time_counter)=spring_length(time_counter,:)*w_out(2:end)+w_out(1);
%     
% end

%%
% dynamic system with memory (3)
 
% f1=2.11;
% f2=3.73;
% f3=4.33;
% f4=5.98;
% 
% 
% y_0=0;
% y=[y_0];
% time_span=0:.001:100;
% u_input=sin(2*pi*f1*time_span).*sin(2*pi*f2*time_span).*sin(2*pi*f3*time_span);
% for i=1:length(time_span)
%     
%    y_t=y(end)+2*(u_input(i)).^3; 
%    y=[y;y_t]; 
%        
% end
% 
% target_function=(.1)*y(2:end);
% 
% Target=target_function(50000:75000);
% 
% 
% w_out=regress(Target,L);
% 
% for time_counter=1:length(T)
%     
%     output_signal(time_counter)=spring_length(time_counter,:)*w_out(2:end)+w_out(1);
%     
% end


%%
% pure sinusoidal function with different input u(t)

% time_vec = 0:0.001:100;
% 
% target_function=transpose(sin(2*pi*time_vec));
% 
% Target=target_function(50000:75000);
% 
% 
% w_out=regress(Target,L);
% 
% for time_counter=1:length(T)
%     
%     output_signal(time_counter)=spring_length(time_counter,:)*w_out(2:end)+w_out(1);
%      
% end

%%
% same target as the input u(t)

time_vec = 0:0.001:100;


f1=2.11;
f2=3.73;
f3=4.33;


amp_target=2.5;
target_function=amp_target*transpose(sin(2*pi*f1*time_vec).*sin(2*pi*f2*time_vec).*sin(2*pi*f3*time_vec));
%target_function=amp_target*transpose(sin(2*pi*f1*time_vec).*sin(2*pi*f2*time_vec).*sin(2*pi*f3*time_vec)+randn(size(time_vec))/50);


Target=target_function(50000:75000);


w_out=regress(Target,L);

for time_counter=1:length(T)
    
    output_signal(time_counter)=spring_length(time_counter,:)*w_out(2:end)+w_out(1);
    
end

%%


figure
plot(T(50000:75000),output_signal(50000:75000),'b--','LineWidth',2)
hold on

plot(T(75000:end),output_signal(75000:end),'k--','LineWidth',2)
hold on


plot(T,target_function,'r','LineWidth',1)
xlim([50,100])







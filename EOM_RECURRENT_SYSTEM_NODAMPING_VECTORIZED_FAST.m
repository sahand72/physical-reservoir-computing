function dpdt = EOM_RECURRENT_SYSTEM_NODAMPING_VECTORIZED_FAST(t,p,all_springs,input_nodes,num_mass,left_mass_index,right_mass_index,amp_target)
% Building the system of equations of motion
% mass index i
% input nodes contains the weights associated to each node. zero for the
% non-input nodes.

f1=2.11;
f2=3.73;
f3=4.33;

%u=amp_target*sin(2*pi*f1*t)*sin(2*pi*f2*t)*sin(2*pi*f3*t); % For other general tasks
%u=amp_target*(sin(2*pi*f1*t)+randn(size(t))/50); % For LP Emulation Task
u=amp_target*(sin(2*pi*f1*t)*sin(2*pi*f2*t)*sin(2*pi*f3*t)); % For emulation task with noise

dpdt = zeros(4*num_mass,1);

for i=1:num_mass

    if i==left_mass_index || i==right_mass_index
       
            dpdt(4*i-3,1) = 0;
            dpdt(4*i-2,1) = 0;
            dpdt(4*i-1,1) = 0;
            dpdt(4*i,1) = 0;
        
    else
    
    mass_i_connections=find(all_springs(:,1)==i);
    connected_springs=all_springs(mass_i_connections,2);
    stiffness_connected=all_springs(mass_i_connections,4:7);
    l0_connections=all_springs(mass_i_connections,3);
    
    w_in=input_nodes(i,2);
              
        
        k1=stiffness_connected(:,1);
        k3=stiffness_connected(:,2);
        %d1=stiffness_connected(:,3);
        %d3=stiffness_connected(:,4);
        
        
        x_connected=sqrt((p(4*i-3)-p(4*connected_springs(:)-3)).^2+(p(4*i-1)-p(4*connected_springs(:)-1)).^2)-l0_connections;

        
        p_x1=-(k3(:).*x_connected(:).^3+k1(:).*x_connected(:)); 
          
        p_x1_vec=p_x1(:).*([p(4*i-3)-p(4*connected_springs(:)-3),p(4*i-1)-p(4*connected_springs(:)-1)]./norm([p(4*i-3)-p(4*connected_springs(:)-3),p(4*i-1)-p(4*connected_springs(:)-1)]));
        
        
        x_spring_force=p_x1_vec(:,1);
        y_spring_force=p_x1_vec(:,2);
        
        
        for j=1:length(connected_springs)
            
        d1=stiffness_connected(j,3);
        d3=stiffness_connected(j,4);
            
        relative_vel(j,:)=[p(4*i-2),p(4*i)]-[p(4*connected_springs(j)-2),p(4*connected_springs(j))];
        
        q_x2_vec=(-d1*(dot(relative_vel(j,:),[p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)])/...
            norm([p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)]))-...
            d3*(dot(relative_vel(j,:),[p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)])/...
            norm([p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)]))^3)*...
            ([p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)]/...
            norm([p(4*i-3)-p(4*connected_springs(j)-3),p(4*i-1)-p(4*connected_springs(j)-1)]));
        x_damping_force(j)=q_x2_vec(1);
        y_damping_force(j)=q_x2_vec(2);
   
        end
        

    
    F_x=sum(x_spring_force)+sum(x_damping_force);
    F_y=sum(y_spring_force)+sum(y_damping_force);
    


    m=1;
    
    dpdt(4*i-3,1) = p(4*i-2);
    dpdt(4*i-2,1) = (1/m)*(F_x+w_in*u);
    dpdt(4*i-1,1) = p(4*i);
    dpdt(4*i,1) = (1/m)*(F_y);
    
    end
    

end
end
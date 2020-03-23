%% Simple toy model

% Simulate the effects of social distancing and what happens when its
% lifted

% This code simulates the transmission of a virus (e.g., COVID-19) due to
% the interaction of people represented as colliding circles. The code
% outputs a gif to view result later. 

% See another type of toy model that simulates COVID transmission for
% people going to the market every so often: github.com/seismotologist & @seismotologist

% Please read the excellent article in the Washting Post by Harry Stevens (washingtonpost.com/graphics/2020/world/corona-simulator/)
% that explains the effects of social distancing which the code was inspired from .

% ------------------------------------------------------
% Author: Chris Milliner, JPL, Caltech, 
% christopher.milliner@jpl.nasa.gov
% https://science.jpl.nasa.gov/people/Milliner/

% Log: Created 03/22/20
% ------------------------------------------------------


clc
clear
% close all

%% Parameters
Num_people              = 200;          % total number of people
Perc_start_sick         = 2;            % Percent of people who start off as sick
Perc_social_distance    = 110;           % percent of people who are using social distancing (i.e., they dont move)
Velocity_factor         = 0.3;          % factor to change velocity of motion
Duration_of_recovery    = 90;          % amount of time it takes to reach recovery in units of number of time
ball_radius             = 0.02;         % ball width, affects liklihood of interacting with another, wider the ball, more likely it will collide with another
Num_time_steps          = 300;          % total number of time steps 
dt                      = 1;            % units of time. Can improve fidelity of simulation. Avoid large values, as this could lead to missing collisions
Time_stop_Social_dist   = 0.5;          % Relative time (between 0 and 1) when to stop social distancing. Set to 1 if you don't want to enact this. 

% Domain extent
left_extent_domain      = -3;           % left extent...
right_extent_domain     = 3;
top_extent_domain       = 1.5;
bot_extent_domain       = -1.5;

% Save result to a gif
Write_output_2_file     = 1;            % switch to turn on (1 = will write out a gif of results).
GIF_delay_time          = 0;            % add a delay time between each frame 
COVID_vid_filename      = './GIF_output/COVID_sim_1.gif';% filename of the gif

% for plotting
Mksize          = 6;                    % size of dots, only for plotting. not for interaction
%---------------------------------------------------------

domain_width = abs(left_extent_domain - right_extent_domain);
domain_height = abs(top_extent_domain - bot_extent_domain);

% Randomize positions
x_vecs = rand(Num_people,1)*domain_width-domain_width/2;
y_vecs = rand(Num_people,1)*domain_height-domain_height/2;

% randomize velocities with uniform distribution
velo_vecs_x = rand(Num_people,1)*Velocity_factor;
velo_vecs_y = rand(Num_people,1)*Velocity_factor;

%-------------------
% Social distancing 
%-------------------
rand_num_peeps_2_SD = floor(length(velo_vecs_x)*(Perc_social_distance/100));% randomize who is using s.d.
rand_indx_peeps_SD = randi([1 length(velo_vecs_x)],rand_num_peeps_2_SD,1);%(TO IMPROVE) this doesnt create a unique set of indices!

% not going anywhere...
velo_vecs_x(rand_indx_peeps_SD) = 0;% people with s.d. stay put
velo_vecs_y(rand_indx_peeps_SD) = 0;

%% Prep vars
Healthy = [ones(Num_people,1),zeros(Num_people,1),zeros(Num_people,1)];% variable to track who is sick, healthy or recovered

rand_peeps_that_r_sick = ceil(length(velo_vecs_x)*(Perc_start_sick/100));% randomize those who are sick
rand_indx_peeps_r_sick = randi([1 length(velo_vecs_x)],rand_peeps_that_r_sick,1);%(TO IMPROVE) this doesnt create a unique set of indices!

Health_tracker(rand_indx_peeps_r_sick) = true;% make them sick

Duration_sick = zeros(Num_people,1);
Duration_sick(rand_indx_peeps_r_sick) = true;
Person_recovered = false(Num_people,1);

Fig_ID          = ceil(rand(1)*1000);       % figure number
Currently_Sick  = zeros(1,length(1:dt:Num_time_steps));
Total_Recovered = zeros(1,length(1:dt:Num_time_steps));
Unaffected      = zeros(1,length(1:dt:Num_time_steps));
k_counter       = 1;

if ~exist('GIF_output', 'dir')
    mkdir GIF_output
end

First_time = 0;
for time_t = 1:dt:Num_time_steps
     
    % update positions
    x_vecs = x_vecs + velo_vecs_x*dt;
    y_vecs = y_vecs + velo_vecs_y*dt;
    
    %----------------------
    % Boundaries
    %----------------------
    
    % Adjust Velocities 
    velo_vecs_y(y_vecs>=top_extent_domain) = velo_vecs_y(y_vecs>=top_extent_domain)*-1;% hit top boundary, go back adjust position + velocity vector     
    velo_vecs_y(y_vecs<=bot_extent_domain) = velo_vecs_y(y_vecs<=bot_extent_domain)*-1;% hit bott. boundary 
    velo_vecs_x(x_vecs<=left_extent_domain) = velo_vecs_x(x_vecs<=left_extent_domain)*-1;% hit left
    velo_vecs_x(x_vecs>=right_extent_domain) = velo_vecs_x(x_vecs>=right_extent_domain)*-1;% hit right 
   
    
    %----------------------
    % Collision detection 
    %----------------------
    dX_matrix = repmat(x_vecs,1,length(x_vecs)) - repmat(x_vecs,1,length(x_vecs))';
    dY_matrix = repmat(y_vecs,1,length(y_vecs)) - repmat(y_vecs,1,length(y_vecs))';
    
    Distances_matrix_between_each_ball = sqrt(dX_matrix.^2 + dY_matrix.^2);% calcualte distances between all particles
    [Indx_collided, ball_2] = find(Distances_matrix_between_each_ball > 0 & Distances_matrix_between_each_ball < 2*ball_radius);% find all balls where the distances is less than twice the radius, therefore a collision has occurred. 
    % note this will miss collisions that may have occurred between the
    % time step. Therefore best to avoid using a dt that is too large.
    
    % Update health
    Indx_r_sick = find(Health_tracker);
    indx_met_other_sick = ismember(ball_2,Indx_r_sick);% find all the collisions where the other person was sick
    Indx_collied_w_sick = Indx_collided(indx_met_other_sick);
    Health_tracker(Indx_collied_w_sick) = true;% make someone sick if they collided when the other person was sick
    Health_tracker(Person_recovered) = false;% turn off if a person had already recovered cannot get sick again. 

    Duration_sick(Health_tracker) = Duration_sick(Health_tracker)+dt;% update time person has been sick
    Person_recovered_now = Duration_sick>=Duration_of_recovery;% person has now recovered
    Person_recovered(Person_recovered_now) = true;% now permanentaly add the recovered person to the recovered list
    
    
    %----------------------
    % Static bumping 
    %----------------------
    if ~isempty(Indx_collided)
        
        Total_dist_between_collided_balls = diag(Distances_matrix_between_each_ball(Indx_collided, ball_2));% amount of overlap (D)
        D_between_collided_balls = 2*ball_radius - Total_dist_between_collided_balls;% Dist of overlap
        Push_back = D_between_collided_balls/2;% scalar amount to push balls apart from each other (D/2)
        
        unit_vector_normal_ball_x = (x_vecs(Indx_collided) - x_vecs(ball_2))./Total_dist_between_collided_balls;% normal vector of direction to push them away from each other.
        unit_vector_normal_ball_y = (y_vecs(Indx_collided) - y_vecs(ball_2))./Total_dist_between_collided_balls;
        
        x_vecs(Indx_collided) = x_vecs(Indx_collided) + unit_vector_normal_ball_x.*Push_back;% update positions
        y_vecs(Indx_collided) = y_vecs(Indx_collided) + unit_vector_normal_ball_y.*Push_back;
    
        unit_vector_normal_ball_x = (x_vecs(Indx_collided) - x_vecs(ball_2))./(2*ball_radius);% update 
        unit_vector_normal_ball_y = (y_vecs(Indx_collided) - y_vecs(ball_2))./(2*ball_radius);
        
        %----------------------
        % Dynamic bumping
        %----------------------
       
        % tangential vector
        unit_vector_tang_ball_x = unit_vector_normal_ball_y*-1;
        unit_vector_tang_ball_y = unit_vector_normal_ball_x;
        
        veloc_in_tang_dir = velo_vecs_x(Indx_collided).*unit_vector_tang_ball_x + velo_vecs_y(Indx_collided).*unit_vector_tang_ball_y;% dot of tangent unit and velocity vector. how much of velocity
        % goes into tangent direction?

        velo_in_normal_dir = zeros(Num_people,1);
        velo_in_normal_dir(Indx_collided) = velo_vecs_x(Indx_collided).*unit_vector_normal_ball_x + velo_vecs_y(Indx_collided).*unit_vector_normal_ball_y; % dot of normal unit and velocity vector. how much of velocity
        % goes into normal direction?

        Indx_is_SD = ismember(ball_2,rand_indx_peeps_SD);% find any other ball indx that is also an SD index.
        % sep those indices that are and are not SD. use below
        other_ball_is_SD = Indx_collided(Indx_is_SD);% this ball is a SD. 
        ball_not_SD = ball_2(~Indx_is_SD);% this ball is not SD
        
        velo_in_normal_dir(ball_not_SD) = velo_in_normal_dir(Indx_collided(~Indx_is_SD));% only applies if other ball is not in SD
        velo_in_normal_dir(other_ball_is_SD) = velo_in_normal_dir(other_ball_is_SD);% if other ball is in SD, then one that isnt should treat other ball as boundary and bounce off.  
        
        velo_vecs_x(Indx_collided) = unit_vector_tang_ball_x.*veloc_in_tang_dir + unit_vector_normal_ball_x.*velo_in_normal_dir(Indx_collided);% New velocity
        velo_vecs_y(Indx_collided) = unit_vector_tang_ball_y.*veloc_in_tang_dir + unit_vector_normal_ball_y.*velo_in_normal_dir(Indx_collided);

        
    end
    
    if time_t/Num_time_steps < Time_stop_Social_dist+1e-12
        % people in SD stay there.
        velo_vecs_x(rand_indx_peeps_SD) = 0;
        velo_vecs_y(rand_indx_peeps_SD) = 0;
        
    elseif time_t/Num_time_steps >= Time_stop_Social_dist && ~First_time
        % first time, then need to randomize velocities of those just
        % realeased from S.D.
        velo_vecs_x(rand_indx_peeps_SD) = rand(length(rand_indx_peeps_SD),1)*Velocity_factor;
        velo_vecs_y(rand_indx_peeps_SD) = rand(length(rand_indx_peeps_SD),1)*Velocity_factor;
        
        First_time = 1;       
        rand_indx_peeps_SD = [];% remove vars so wont affect collision 
        
    end
    
    % Count up 
    Currently_Sick(k_counter) = sum(double(Health_tracker));% total sick
    Total_Recovered(k_counter) = sum(double(Person_recovered));% total recovered at the moment.
    Unaffected(k_counter) = Num_people - Currently_Sick(k_counter) - Total_Recovered(k_counter) + Currently_Sick(k_counter);% current total number unaffected (neither of the above), add the sick in again for plotting
    
    
    %----------------------
    % Plotting
    %----------------------
    
    h = figure(Fig_ID);
    clf(Fig_ID)
    subplot(8,4,[2,3,6,7])
    
    % plot recovered
    fill([0,1:dt:time_t,time_t],[0,ones(size(Unaffected(1:k_counter)))*Num_people,0],[109, 140, 252]/256)
    hold on 
    plot([1:dt:time_t],ones(size(Unaffected(1:k_counter)))*Num_people,'Color',[109, 140, 252]/256,'linewidth',2)
    % plot unaffected
    fill([0,1:dt:time_t,time_t],[0,Unaffected(1:k_counter),0],[47, 255, 189]/256)
    plot([1:dt:time_t],Unaffected(1:k_counter),'Color',[47, 255, 189]/256,'linewidth',2)
    % plot sick 
    fill([0,1:dt:time_t,time_t],[0,Currently_Sick(1:k_counter),0],[235, 102, 0]/256)
    plot([1:dt:time_t],Currently_Sick(1:k_counter),'Color',[235, 102, 0]/256,'linewidth',2)

    axis([1 Num_time_steps 0 Num_people])
    set(gca,'fontsize',12)
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
 
    % text
    txt_topleft = 'Healthy';
    txthandle = text(-160,Num_people*0.7,txt_topleft);%,'Color',[47, 255, 189]/256);
    set(txthandle,'fontsize',12)%,'FontWeight','bold')
    
    txt_middle = 'Sick';
    txthandle = text(-160,Num_people*0.5,txt_middle);%,'Color',[235, 102, 0]/256);
    set(txthandle,'fontsize',12)%,'FontWeight','bold')
    
    txt_topleft = 'Recovered';
    txthandle = text(-160,Num_people*0.3,txt_topleft);%,'Color',[109, 140, 252]/256);
    set(txthandle,'fontsize',12)%,'FontWeight','bold')
    
    % data
    txt_topleft = num2str(Unaffected(k_counter) - Currently_Sick(k_counter));
    txthandle = text(-40,Num_people*0.7,txt_topleft,'Color',[47, 255, 189]/256);
    set(txthandle,'fontsize',12,'FontWeight','bold')
    
    txt_middle = num2str(Currently_Sick(k_counter));
    txthandle = text(-40,Num_people*0.5,txt_middle,'Color',[235, 102, 0]/256);
    set(txthandle,'fontsize',12,'FontWeight','bold')
    
    txt_topleft = num2str(Total_Recovered(k_counter));
    txthandle = text(-40,Num_people*0.3,txt_topleft,'Color',[109, 140, 252]/256);
    set(txthandle,'fontsize',12,'FontWeight','bold')
    
    
    
    subplot(8,4,[9:32])
    % plot all (will be plotted as unaffected)
    plot(x_vecs,y_vecs,'o','markerfacecolor',[47, 255, 189]/256,'markeredgecolor',[47, 255, 189]/256,'markersize',Mksize)
    hold on 
    % plot sick
    plot(x_vecs(Health_tracker),y_vecs(Health_tracker),'o','markerfacecolor',[235, 102, 0]/256,'markeredgecolor',[235, 102, 0]/256,'markersize',Mksize)
    % plot recovered
    plot(x_vecs(Person_recovered),y_vecs(Person_recovered),'o','markerfacecolor',[109, 140, 252]/256,'markeredgecolor',[109, 140, 252]/256,'markersize',Mksize)
    axis equal
    axis([left_extent_domain right_extent_domain bot_extent_domain top_extent_domain])
    set(gca,'fontsize',15)
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gcf,'color','w');
    set(gca,'linewidth',2)% figure boundary 
    
    drawnow 
    switch Write_output_2_file
        case 1
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            
            if time_t == 1
                imwrite(imind,cm,COVID_vid_filename,'gif', 'Loopcount',inf,'DelayTime',GIF_delay_time);
            else
                imwrite(imind,cm,COVID_vid_filename,'gif','WriteMode','append','DelayTime',GIF_delay_time);
            end
    end
    
    k_counter = k_counter +1;
end

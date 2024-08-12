classdef steady_diffusionModel_toy
    %DIFFUIONMODEL Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods


        % return 1 if a pos_x pos_y fallis inside polygon defined by bounding_x bounding_y
        function is_inside = is_point_inside_polygon(obj , pos_x , pos_y, bounding_x , bounding_y)

            num_vertices = size(bounding_x, 2);

            % Initialize the crossing count
            crossing_count = 0;

            % Iterate through each edge of the polygon
            for i = 1:num_vertices
                x1 = bounding_x(i);
                y1 = bounding_y(i);
                if i == num_vertices
                    x2 = bounding_x(1);
                    y2 = bounding_y(1);
                else
                    x2 = bounding_x(i + 1);
                    y2 = bounding_y(i + 1);
                end

                % Check for intersection between the horizontal ray and the edge
                if (pos_y > min(y1, y2)) && (pos_y <= max(y1, y2)) && (pos_x <= max(x1, x2)) && (y1 ~= y2)
                    % Calculate the x-coordinate of the intersection point
                    x_intersection = (pos_y - y1) * (x2 - x1) / (y2 - y1) + x1;

                    % Check if the intersection point is to the right of the test point
                    if x_intersection > pos_x
                        crossing_count = crossing_count + 1;
                    end
                end
            end

            % If the crossing count is odd, the point is inside the polygon
            is_inside = mod(crossing_count, 2) == 1;
        end

        % find the point of intersection and angle of incidence for two vectors; return collision = 0 if the vectors don't intersect
        function  [collision , i_x , i_y , angle_in] = get_line_intersection(obj, p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y)
            i_x = 0;
            i_y = 0;
            angle_in = 0;
            s1_x = p1_x - p0_x;
            s1_y = p1_y - p0_y;
            s2_x = p3_x - p2_x;
            s2_y = p3_y - p2_y;

            s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
            t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

            if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
                % Collision detected
                i_x = p0_x + (t * s1_x);
                i_y = p0_y + (t * s1_y);
                collision = 1;
                dp = dot([(p1_x - p0_x) , (p1_y - p0_y)] , [(p3_x - p2_x) , (p3_y - p2_y)]);
                angle_in = acos(dp  / (norm([(p1_x - p0_x) , (p1_y - p0_y)]) * norm([(p3_x - p2_x) , (p3_y - p2_y)]) ));

            else
                % No collision
                collision = 0;
            end
        end

        % update a particle's position with reflective bouncing at polygon boundaires defined by bounding_x bounding_y
        function [new_pos_x , new_pos_y , new_delta_x , new_delta_y , intersect_side] = update_position_reflective(obj , old_x , old_y , bounding_x , bounding_y , delta_x , delta_y , prev_intersect , do_plot)

        %disp(["size of bounding_x should be n_vertices: " size(bounding_x,2)])
            intersect_side = 0;
            % for every line segement, check for intersection
            found_intersection = false; % helper
            for side = 1:size(bounding_x,2)
                if found_intersection == false && side ~= prev_intersect
                    this_point = side;
                    next_point = side+1;
                    % wrap around point list
                    if next_point > size(bounding_x,2)
                        next_point = 1;
                    end
                        
                    if do_plot == true
                        plot([old_x , old_x+delta_x] , [old_y , old_y+delta_y], 'Color','y' , 'LineWidth' , 2) % plot the  target displacement
                    end
                    % check for collision
                    [col , ix , iy , ang] = obj.get_line_intersection(old_x,old_y , old_x+delta_x,old_y+delta_y , bounding_x(this_point),bounding_y(this_point), bounding_x(next_point),bounding_y(next_point));

                    % if it collides, calculate bounce
                    if col == 1
                        found_intersection = true;
                       % disp('Found intersection set to true')
                        intersect_side = side;
                        if do_plot == true
                            plot([old_x , ix] , [old_y , iy], 'Color','r' , 'LineWidth'  , 2) % plot from the start to the intersection
                        end

                        % this is janky and could likely be cleaned up
                        %calculate new starting point, direction & magnitude left to travel
                        mag_traveled = norm([ix-old_x , iy-old_y]);
                        mag_remaining = norm([delta_x , delta_y]) - mag_traveled;

                        edge_delta_x = bounding_x(next_point) - bounding_x(this_point);
                        edge_delta_y = bounding_y(next_point) - bounding_y(this_point);
                        edge_normal = [-edge_delta_y edge_delta_x];
                        edge_normal = edge_normal / norm(edge_normal);

                        reflected_vector = [delta_x delta_y] - 2 * ([delta_x delta_y] * (edge_normal')) * edge_normal;

                        reflected_unit_vector = reflected_vector / norm(reflected_vector);

                        % calculate new starting point at intersection point &
                        % new delta's based on angle of incidence
                        new_pos_x = ix;
                        new_pos_y = iy;

                        new_delta_x = reflected_unit_vector(1)*mag_remaining;
                        new_delta_y = reflected_unit_vector(2)*mag_remaining;
                        
                        %disp(["reflection angle is " , reflection_angle]);
                        %disp(["Side hit is " , side]);
                    end
                end
            end
            % if it checks every side and doesn't flag, update pos_x
            if found_intersection == false
                new_pos_x = old_x + delta_x;
                new_pos_y = old_y + delta_y;
                new_delta_x = 0;
                new_delta_y = 0;
            end

        end

        % wrapper function to call update_position_reflective until there is no more reflection
        function [pos_x , pos_y] = wrap_update_position_reflective(obj , pos_x , pos_y , delta_x , delta_y , bounding_x , bounding_y , do_plot)
            intersect_side = 0;
            while delta_x ~= 0 && delta_y ~= 0
                %plot([pos_x,pos_x+delta_x] ,[pos_y,pos_y+delta_y] , 'Color' ,'r' , LineWidth=2)
                [pos_x , pos_y , delta_x , delta_y ,intersect_side] = obj.update_position_reflective(pos_x , pos_y , bounding_x , bounding_y , delta_x , delta_y , intersect_side , do_plot);
            end
        end


        function [x_vertices , y_vertices] = generate_polygon_field(obj , square_radius , gridSize , n_vertices)

            % Create figure
            %figure;
            %axis equal;
            %hold on;

            x_vertices = zeros(gridSize * gridSize , n_vertices);
            y_vertices = zeros(gridSize* gridSize , n_vertices);
            angle_intervals = (360/n_vertices/2):(360/n_vertices):360;

            proto_shape_x = [0 + square_radius * cosd(angle_intervals)];
            proto_shape_y = [0 + square_radius * sind(angle_intervals)];

            width = max(proto_shape_x) - min(proto_shape_x);
            height = max(proto_shape_y) - min(proto_shape_y);
            side_length = proto_shape_y(1)*2;

            % Loop through rows and columns to create hexagons
            for i = 1:gridSize
                y_offset = (i-1) * side_length / 2;
                for j = 1:gridSize
                    x_offset = (j-1) * side_length / 2;
                    this_poly = ((i-1)*gridSize+j);
                    x_vertices(this_poly,:) = x_offset + proto_shape_x;
                    y_vertices(this_poly,:) = y_offset + proto_shape_y;

                    %fill([x_vertices(this_poly,:)] , y_vertices(this_poly,:) , 'b' , 'EdgeColor' , 'k')

                end
            end

            %hold off;
        end

        % return a new particle with position randomly assigned within polygon defined by bounding_x and bounding_y vertices, and popoulation drawn from initial_probs
        function [x y pop] = initialize_agent(obj , bounding_x , bounding_y , initial_probs)
            inside_target_cell = 0;

            while inside_target_cell == 0
                a = min(bounding_x);
                b = max(bounding_x);
                x = (b-a) * rand(1,1) + a;

                a = min(bounding_y);
                b = max(bounding_y);
                y = (b-a).*rand(1,1) + a;

                inside_target_cell = obj.is_point_inside_polygon(x , y , bounding_x , bounding_y);

            end
            pop = randsample(1:size(initial_probs,2), 1, true, initial_probs);

        end

        function [msds , final_dist_from_sender , fluxxed] = run_simulation(obj , n_agents , max_time , grid_x , grid_y ,rates, transition_matrix , initial_probs , constrained)

            tic;

            pos_x = zeros(n_agents,1);
            pos_y = zeros(n_agents,1);
            pops = zeros(n_agents,1);
            cell = zeros(n_agents,1);

            msds = zeros(max_time,1);
            fluxxed = zeros(max_time,1);
            start_cell = floor((size(grid_x,1)) / 2) + 1;

            home_x = mean([min(grid_x(start_cell,:)) , max(grid_x(start_cell,:))]);
            home_y = mean([min(grid_y(start_cell,:)) , max(grid_y(start_cell,:))]);

            for z = 1:max_time
                if mod(z , 5000) == 0
                    disp(['Iteration: ', num2str(z)])
                end

                % top off missing particles / only runs first iteration unless cells signal/are degraded
                ind_missing_particles = find(pops == 0);
                for i = 1:length(ind_missing_particles)
                    q = ind_missing_particles(i);
                    [pos_x(q) , pos_y(q) , pops(q)] = obj.initialize_agent(grid_x(start_cell,:) , grid_y(start_cell,:) , initial_probs);
                    cell(q) = start_cell;
                end

                for i = 1:n_agents
                    [pos_x(i) , pos_y(i) , cell(i)] = obj.generate_new_coordinates(pos_x(i) , pos_y(i) , pops(i) , cell(i) , grid_x , grid_y , rates , constrained);
                    pops(i) = obj.update_pop(pops(i) , transition_matrix);
                end
                %pops(pops == 6) = 0; % mark particles for removal if they signaled

                % calculate the MSD from the origin
                msds(z) = mean( ((pos_x - home_x).^2 + (pos_y - home_y).^2) );
                fluxxed(z) = length(find(cell ~= start_cell));

            end

            % % for each particle, calculate the distance to cell bounadry
            % % only do this after model hits steady state, to plot the final distribution
            % final_dist_from_sender = zeros(n_agents , 1);
            % for i = 1:n_agents
            %     dist_x = min(abs(pos_x(i) - cell_bounds(start_cell,[1 2]))); % min distance to x boundary of sender cell
            %     dist_y = min(abs(pos_y(i) - cell_bounds(start_cell,[3 4]))); % min distance to y boundary of sender cell
            %     final_dist_from_sender(i) = sqrt(dist_x^2 + dist_y^2);
            % end
            % final_dist_from_sender(cell == start_cell) = -final_dist_from_sender(cell == start_cell); % remove particles that never left home
            final_dist_from_sender = 0;
            toc;
        end




        %%
        %%
        %%
        %%
        % based on anders's problem set 20.430; 2
        % simulate diffusion in 2D by independently simulating
        % diffusion in 2 different dimensions


        function [new_x , new_y , new_cell] = generate_new_coordinates(obj ,old_x , old_y , pop , old_cell , grid_x , grid_y  , rates , constrained)

            signaled = 0;
            %rates = [0 0.2 1 12];

            new_cell = 0;

            switch pop
                case 1 %% pop 1 is either constrained or free
                    D = rates(1);
                    if constrained
                        if old_cell == 0
                            old_cell = obj.update_cell(old_x , old_y , grid_x , grid_y); % for particles thare are *newly* confined give them a cell
                        end
                        tau = 0.006; % 6 ms
                        L = sqrt(2*D*tau);
                        delta_x = normrnd(0,L);
                        delta_y = normrnd(0,L);
                        [new_x , new_y] = obj.wrap_update_position_reflective(old_x , old_y , delta_x , delta_y , grid_x(old_cell,:) , grid_y(old_cell,:) , false);
                        new_cell = old_cell;
                    else
                        tau = 0.006;
                        L = sqrt(2*D*tau);
                        new_x = old_x + normrnd(0,L);
                        new_y = old_y + normrnd(0,L);
                        %new_cell = obj.update_cell(new_x , new_y , grid_x , grid_y);
                        new_cell = 0;
                    end

                case 2 %% pop 2 is always free
                    D = rates(2);
                    tau = 0.006;
                    L = sqrt(2*D*tau);
                    new_x = old_x + normrnd(0,L);
                    new_y = old_y + normrnd(0,L);
                    new_cell = 0;

                otherwise
                    error('disallowed population value')
                    disp('error')
            end

        end

      

        function new_pop = update_pop(obj ,old_pop , transition_matrix)

            transition_vector = transition_matrix(old_pop,:);
            new_pop = randsample(1:size(transition_matrix,1), 1, true, transition_vector);

        end

        function [starting_probs] = generate_intial_probs(obj ,frac_fast)

            % below is based on : pickle_dir = '/Volumes/gschliss/software/spt_tracking_code/pkl_files/2023_03_21/dt=6ms_ld=1p71_q=3p0_windowLength=5_nbSubstep=1_prior=new1p72_splitSize=10_maxLen=3-300_singleCheck'
            starting_probs = [1-frac_fast frac_fast];

        end

        function [matrix] = generate_transition_matrix(obj ,tx_rate , fraction_fast)

            % below is based on : pickle_dir = '/Volumes/gschliss/software/spt_tracking_code/pkl_files/2023_03_21/dt=6ms_ld=1p71_q=3p0_windowLength=5_nbSubstep=1_prior=new1p72_splitSize=10_maxLen=3-300_singleCheck'

            matrix = [0 tx_rate * fraction_fast
                tx_rate 0];

            for i = 1:size(matrix,1)
                for j = 1:size(matrix,2)
                    if(i == j)
                        matrix(i,j) = 1-sum(matrix(i,:));
                    end
                end
            end
        end



        function new_cell = update_cell(obj ,x , y , grid_x , grid_y)
            
            new_cell_options = 0;
            for j = 1:size(grid_x,1)
                if obj.is_point_inside_polygon(x , y , grid_x(j,:) , grid_y(j,:))
                    new_cell_options = [new_cell_options , j];
                end
            end
            if new_cell_options == 0
                disp('cell finding error')
            end
            new_cell_options = new_cell_options(new_cell_options>0);
            new_cell = randsample(new_cell_options , 1);

        end


    end
end


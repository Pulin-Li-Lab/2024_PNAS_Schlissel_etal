classdef steady_diffusionModel_toy
    %DIFFUIONMODEL Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods


        function [x y pop] = initialize_agent(obj , cell_bounds , start_cell , initial_probs)
            a = min(cell_bounds(start_cell,1:2));
            b = max(cell_bounds(start_cell,1:2));
            x = (b-a).*rand(1,1) + a;

            a = min(cell_bounds(start_cell,3:4));
            b = max(cell_bounds(start_cell,3:4));
            y = (b-a).*rand(1,1) + a;

            pop = randsample(1:size(initial_probs,2), 1, true, initial_probs);

        end

        function [pos_x , pos_y , pops , cell] = initialize_vectors(obj , n_agents);
            pos_x = zeros(n_agents,1);
            pos_y = zeros(n_agents,1);
            pops = zeros(n_agents,1);
            cell = zeros(n_agents,1);
        end


        function [msds , final_dist_from_sender , fluxxed , backwards] = run_simulation(obj , sender_conc , max_time , cell_bounds ,rates, transition_matrix , initial_probs , constrained)
            tic;

            [send_pos_x , send_pos_y, send_pops , send_cell] = obj.initialize_vectors(sender_conc);
            send_pops = randsample([1 2],length(send_pops), true , initial_probs)';

            % initialize empty vectors for molecules in the receiver field
            [pos_x , pos_y , pops , cell] = obj.initialize_vectors(0);

            msds = zeros(max_time,1); % initialize empty msds output vector for MSD at each time step
            fluxxed = zeros(max_time, 1); %% initialize vector to store number of molecules secreted over time
            start_cell = floor((size(cell_bounds,1)) / 2) + 1; % identify start cell in the middle of the field

            backwards = zeros(max_time, 1);

            home_x = mean([min(cell_bounds(start_cell,1:2)) , max(cell_bounds(start_cell,1:2))]); %% identify center of start cell
            home_y = mean([min(cell_bounds(start_cell,3:4)) , max(cell_bounds(start_cell,3:4))]);

            for z = 1:max_time

                % report progress
                if mod(z , 1000) == 0
                    disp(['Iteration: ', num2str(z)])
                end


                % every time step, check to see if the sender population is "full", and if it isn't replace the
                % missing molecules
                ind_missing_sender = find(send_pops == 6 | send_cell ~= start_cell);
                for i = 1:length(ind_missing_sender)
                    q = ind_missing_sender(i);
                    [send_pos_x(q) , send_pos_y(q) , dummy] = obj.initialize_agent(cell_bounds , start_cell , initial_probs); % allow send_pops(q) to re-use the previous value
                    send_cell(q) = start_cell;
                end
        		send_pops = randsample([1 2],length(send_pops), true , initial_probs)'; % force sender cell population distribution to match initial_probs at every time step

                % let all sender molecules move (but don't do anything with them yet)

                for i = 1:length(send_pos_x)
                    [send_pos_x(i) , send_pos_y(i) , send_cell(i)] = obj.generate_new_coordinates(send_pos_x(i) , send_pos_y(i) , send_pops(i) , send_cell(i) , cell_bounds , rates , constrained);
                    %send_pops(i) = obj.update_pop(send_pops(i) , transition_matrix);
                end
                to_transfer = find(send_cell ~= start_cell);


                % let all of the molecules in the field move
                for i = 1:length(pos_x)
                    [pos_x(i) , pos_y(i) , cell(i)] = obj.generate_new_coordinates(pos_x(i) , pos_y(i) , pops(i) , cell(i) , cell_bounds , rates , constrained);
                    pops(i) = obj.update_pop(pops(i) , transition_matrix);
                end


        		% any molecules from sender vector that left the sender cell, transfer to the corresponding main molecule vectors
                pos_x((end+1):(end+length(to_transfer))) = send_pos_x(to_transfer);
                pos_y((end+1):(end+length(to_transfer))) = send_pos_y(to_transfer);
                cell((end+1):(end+length(to_transfer))) = send_cell(to_transfer);
                pops((end+1):(end+length(to_transfer))) = send_pops(to_transfer);

                backwards(z) = length(find(cell == start_cell));
                % delete any molecule that re-entered the start cell
                to_remove = find(cell == start_cell);

                % doing it this way will basically be like re-initailizing molecules near the border where they jumped backwards.
        		% this will fail if the number of molecules jumping backwards is greater than the sender concentration
                for i = 1:length(to_remove)
        		    if i > size(send_pos_x,1)
                        disp('removed too many particles')
                        continue
                    end
                    send_pos_x(i) = pos_x(to_remove(i)); %replace first (i.e. arbitrary) particles with backwards-jumped particle
                    send_pos_y(i) = pos_y(to_remove(i));
                    send_cell(i) = start_cell;
                    send_pops(i) = pops(to_remove(i));
                end

                pos_x(to_remove) = [];
                pos_y(to_remove) = [];
                cell(to_remove) = [];
                pops(to_remove) = [];

                % calculate the MSD from the origin for all of the molecules in the receiver field
                msds(z) = mean(((pos_x - home_x).^2 + (pos_y - home_y).^2));

                % record the number of molecules in the receiver field
                fluxxed(z) = length(pos_x);

            end

            % for each particle, calculate the distance to cell bounadry
            % only do this at the final time step, to plot the final distribution
            final_dist_from_sender = zeros(length(pos_x) , 1);
            for i = 1:length(pos_x)
                dist_x = min(abs(pos_x(i) - cell_bounds(start_cell,[1 2]))); % min distance to x boundary of sender cell
                dist_y = min(abs(pos_y(i) - cell_bounds(start_cell,[3 4]))); % min distance to y boundary of sender cell
                final_dist_from_sender(i) = sqrt(dist_x^2 + dist_y^2);
            end
            final_dist_from_sender(cell == start_cell) = -final_dist_from_sender(cell == start_cell); % remove particles that never left home
            toc;
        end




        %%
        %%
        %%
        %%
        % based on anders's problem set 20.430; 2
        % simulate diffusion in 2D by independently simulating
        % diffusion in 2 different dimensions


        function [new_x , new_y , new_cell] = generate_new_coordinates(obj ,old_x , old_y , pop , old_cell , cell_bounds  , rates , constrained)

            signaled = 0;
            %rates = [0 0.2 1 12];

            switch pop
                case 1 %% pop 1 is either constrained or free
                    D = rates(1);
                    if constrained
                        % constrained generator, D = 0.2 & boundary
                        [new_x , new_y] = obj.constrained_generator(old_x , old_y , D , ...
                            cell_bounds(old_cell,1) , cell_bounds(old_cell,2) , cell_bounds(old_cell,3) , cell_bounds(old_cell,4));
                        new_cell = old_cell;
                    else
                        [new_x, new_y] = obj.unconstrained_generator(old_x , old_y , D);
                        new_cell = obj.update_cell(new_x , new_y , cell_bounds);
                    end

                case 2 %% pop 2 is always free
                    D = rates(2);
                    [new_x, new_y] = obj.unconstrained_generator(old_x , old_y , D);
                    new_cell = obj.update_cell(new_x , new_y , cell_bounds);

                otherwise
                    error('disallowed population value')
            end

        end

        function [new_x , new_y] = unconstrained_generator(obj ,old_x , old_y , D)
            tau = 0.006;
            L = sqrt(2*D*tau);
            new_x = old_x + normrnd(0,L);
            new_y = old_y + normrnd(0,L);
        end

        function [new_x , new_y] = constrained_generator(obj ,old_x , old_y , D , x_min , x_max , y_min , y_max)

            tau = 0.006; % 6 ms

            L = sqrt(2*D*tau);

            target_x = old_x + normrnd(0,L);
            while(target_x < x_min || target_x > x_max)
                if target_x < x_min %% if it falls out left
                    how_far_over = (x_min - target_x);
                    target_x = x_min + how_far_over;
                else %% if it falls out right
                    how_far_over = x_max - target_x;
                    target_x = x_max + how_far_over;
                end

            end
            new_x = target_x;

            target_y = old_y + normrnd(0,L);
            while(target_y < y_min || target_y > y_max)
                if target_y < y_min  %% if it falls out bottom
                    how_far_over = (y_min - target_y);
                    target_y = y_min + how_far_over;
                else %% if it falls out top
                    how_far_over = y_max - target_y;
                    target_y = y_max + how_far_over;
                end
            end
            new_y = target_y;

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



        function new_cell = update_cell(obj ,x , y , cell_bounds)
    	    if x < 0 || x > max(max(cell_bounds)) || y < 0 || y > max(max(cell_bounds))
                disp('Molecule left the field')
            end
            poss_cell = zeros(size(cell_bounds));
            poss_cell(:,1) = cell_bounds(:,1) <= x;
            poss_cell(:,2) = cell_bounds(:,2) > x;
            poss_cell(:,3) = cell_bounds(:,3) <= y;
            poss_cell(:,4) = cell_bounds(:,4) > y;
            [~ , new_cell] = max(sum(poss_cell , 2));

        end


        % Return a matrix with 4 columns & n_cell_x * n_cells_y rows, where
        % each row is the x_min x_max y_min y_max of a cell in a  grid
        % the cells are assumed to be square, with dimension given by
        % cell_dimension
        function cell_bounds = build_cell_field(obj , n_cells_x , n_cells_y , cell_dimension)

            % Compute the total number of cells
            n_cells = n_cells_x * n_cells_y;

            % Initialize a matrix to store the cell boundaries
            cell_bounds = zeros(n_cells, 4);

            % Compute the x and y coordinates of each cell
            [x, y] = meshgrid(1:n_cells_x, 1:n_cells_y);

            % Compute the left, right, bottom, and top boundaries of each cell
            left_bounds = (x - 1) * cell_dimension;
            right_bounds = x * cell_dimension;
            bottom_bounds = (y - 1) * cell_dimension;
            top_bounds = y * cell_dimension;

            % Store the cell boundaries in the output matrix
            cell_bounds(:, 1) = left_bounds(:);
            cell_bounds(:, 2) = right_bounds(:);
            cell_bounds(:, 3) = bottom_bounds(:);
            cell_bounds(:, 4) = top_bounds(:);

        end


    end
end


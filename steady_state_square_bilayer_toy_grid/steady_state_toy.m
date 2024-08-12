dm =  steady_diffusionModel_toy;

rng(pid)

grates = [slowRate fastRate];

% pass these parameters as arguments
n_agents
max_time
frac_fast
cell_sides
map_fold
map_frac
slowRate
fastRate
constraint

if constraint == 0
    conditions = {'free'};
else
    conditions = {'constrained'};
end
conditions


cell_area

cell_dimension = sqrt((2 * cell_area / cell_sides) / sin(2*pi/cell_sides))
[grid_x grid_y] = dm.generate_polygon_field(cell_dimension,55,4);

%%


%% Test section
% clear
% dm =  steady_diffusionModel_toy;
% obj = dm;
% 
% profile on
% 
% n_agents = 100;
% max_time = 1000;
% frac_fast = 0.05;
% cell_area = 25;
% cell_sides = 4;
% map_fold = 50;
% map_frac = 9;
% slowRate = 0.5;
% fastRate = 1;
% constraint = 1;
% tx_rate = 0.081633;
% 
% 
% grates = [slowRate fastRate];
% 
% if constraint == 0
%     conditions = {'free'};
% else
%     conditions = {'constrained'};
% end
% conditions;
% 
% cell_area
% 
% cell_dimension = sqrt((2 * cell_area / cell_sides) / sin(2*pi/cell_sides))
% [grid_x grid_y] = dm.generate_polygon_field(cell_dimension,55,cell_sides);
% 
% start_cell = floor((size(grid_x,1)) / 2) + 1;
% 
% bounding_x = grid_x(start_cell,:);bounding_y = grid_y(start_cell,:);
% figure; hold on; axis equal; fill(bounding_x , bounding_y , 'g')
% 
% %bounding_x = grid_x(start_cell+1,:);bounding_y = grid_y(start_cell+1,:);fill(bounding_x , bounding_y , 'r')
% 
% pos_x = 136; pos_y = 135; delta_x = -17 ; delta_y = -16; prev_intersect = 0; do_plot = true; intersect_side = 0; old_x = pos_x; old_y = pos_y;
% 
% rates = grates(1,:);
% constrained = strcmp(conditions(1) , 'constrained');
% initial_probs = dm.generate_intial_probs(frac_fast);
% transition_matrix = dm.generate_transition_matrix(tx_rate,frac_fast);
% 
% [msds , final_dist , fluxxed] = dm.run_simulation(n_agents , max_time , grid_x , grid_y , rates , transition_matrix , initial_probs , constrained);


%%


% set up map/reduce; break into map_fold groups and 
% slice tx_rates and take the map_frac group of tx_rates
% this will save all of the files separately, and I can reduce in R
		       
interval = 2; % set by number of processors used
tx_rates = linspace(0, 1 , (map_fold * interval / 2));
start = (map_frac - 1) * (interval/2) + 1;
stop = start + (interval/2) - 1;
if stop > size(tx_rates,2)
    stop = size(tx_rates,2);
end
if start > size(tx_rates,2)
    disp('ERROR -- calling tx_rates that dont exist')
    exit;
end
tx_rates = tx_rates(start:stop)
size(tx_rates)


%parpool(interval)

params = struct();
for n = 1:length(tx_rates)
    this_tx = (tx_rates(n));
    for c = 1:length(conditions)
        this_con = conditions(c);
        con = strcat(this_con,'_', num2str(n));
        con = con{1};
        params.(con).name = con;
        params.(con).tx_rate = this_tx;
        params.(con).conditionIndex = conditions(c);
        params.(con).n_agents = n_agents;
        params.(con).max_time = max_time;
        params.(con).frac_fast = frac_fast;
    end
end
%%

fn = fieldnames(params);
ss_data = {};
for n = 1:length(fn)
    ss_data{n} = fn(n);
end

fn = fieldnames(params);
for n = 1:length(fn)
    f = fn{n};
    tx_rate = params.(f).tx_rate;
    rates = grates(1,:);

    constrained = strcmp(params.(f).conditionIndex{1} , 'constrained');

    initial_probs = dm.generate_intial_probs(frac_fast);

    transition_matrix = dm.generate_transition_matrix(tx_rate,frac_fast);

    [msds , final_dist , fluxxed] = dm.run_simulation(n_agents , max_time , grid_x , grid_y , rates , transition_matrix , initial_probs , constrained);
    
    ss_data{n} = struct();

    ss_data{n}.name = f;
    ss_data{n}.msds = msds;
    ss_data{n}.fluxxed = fluxxed;
    ss_data{n}.final_dist_from_sender = final_dist;
    ss_data{n}.frac_fast = frac_fast;
    ss_data{n}.matrix = transition_matrix;
    ss_data{n}.n_agents = n_agents;
    ss_data{n}.max_time = max_time;

end



index = round(mod(posixtime(datetime('now')),1) * 1000000);

outFile = ['steady_state_test_data/estDRate_toy_agents=' int2str(n_agents) '_max_time=' int2str(max_time) '_cellArea=' int2str(cell_area) '_fracFast=' num2str(frac_fast) '_slowRate=' num2str(slowRate) '_fastRate=' num2str(fastRate) '_condition=' char(conditions) '_map' num2str(map_frac) 'of' num2str(map_fold) '_randIndex_' num2str(index) '_output.mat']
  
save(outFile, 'ss_data')



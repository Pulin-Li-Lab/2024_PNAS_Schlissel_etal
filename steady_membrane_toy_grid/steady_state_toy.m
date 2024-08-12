
dm =  steady_diffusionModel_toy;

%grates = [0.5 20];

grates = [slowRate fastRate];

% pass these parameters as arguments
%sender_conc = 5000;
%max_time = 50000;
%cell_dimension = 5;

%c = parcluster;
%c.NumWorkers = 55;
%c.NumThreads = 2
parpool(12)
  
sender_conc
max_time
frac_fast
cell_dimension
map_fold
map_frac

%%
cell_bounds = dm.build_cell_field(71 , 71 , cell_dimension);


%names = {'zero' , 'ten'}
conditions = {'constrained' , 'free'}
tx_rates = 0:0.01:0.5;
tx_rates = 0:0.002:0.1;
tx_rates = 0.1002:0.002:0.2;	

% set up map/reduce; break into map_fold groups and 
% slice tx_rates and take the map_frac group of tx_rates
% this will save all of the files separately, and I can reduce in R
interval = round(size(tx_rates,2) / map_fold);
start = (map_frac - 1) * interval+1;
stop = start + interval - 1;
if stop > size(tx_rates,2)
    stop = size(tx_rates,2);
end
if start > size(tx_rates,2)
    disp('ERROR -- calling tx_rates that dont exist')
    exit;
end
tx_rates = tx_rates(start:stop)
size(tx_rates)

%%
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
        params.(con).sender_conc = sender_conc;
        params.(con).max_time = max_time;
        params.(con).frac_fast = frac_fast;
    end
end
fn = fieldnames(params);
ss_data = {};
for n = 1:length(fn)
    ss_data{n} = fn(n);
end

fn = fieldnames(params);
parfor n = 1:length(fn)
    f = fn{n};
    tx_rate = params.(f).tx_rate;
    rates = grates(1,:);

    constrained = strcmp(params.(f).conditionIndex{1} , 'constrained');

    initial_probs = dm.generate_intial_probs(frac_fast);

    transition_matrix = dm.generate_transition_matrix(tx_rate,frac_fast);

    [msds , final_dist_from_sender , fluxxed] = dm.run_simulation(sender_conc, max_time , cell_bounds , rates , transition_matrix , initial_probs , constrained);

    ss_data{n} = struct();

    ss_data{n}.name = f;
    ss_data{n}.msds = msds;
    ss_data{n}.fluxxed = fluxxed;
    ss_data{n}.final_dist_from_sender = final_dist_from_sender;
    ss_data{n}.frac_fast = frac_fast;
    ss_data{n}.matrix = transition_matrix;
    ss_data{n}.sender_conc = sender_conc;
    ss_data{n}.max_time = max_time;

end


index = round(mod(posixtime(datetime('now')),1) * 1000000);

outFile = ['steady_state_test_data/estFluxRate_toy_agents=' int2str(sender_conc) '_max_time=' int2str(max_time) '_cellSize=' int2str(cell_dimension) '_fracFast=' num2str(frac_fast) '_slowRate=' num2str(slowRate) '_fastRate=' num2str(fastRate) '_map' num2str(map_frac) 'of' num2str(map_fold) '_randIndex_' num2str(index) '_fixedBackflow_output.mat']

save(outFile, 'ss_data')

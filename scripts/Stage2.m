%===============
% STAGE 2
%===============
if isfile('../aux/log_stage2.csv')
    disp('Stage 2 has already been run! To run once again, delete the file `log_stage2.csv`. Back up old data before rerunning.');
    return;
end
%-------------------
%BAND-RELATED INPUTS
%-------------------
%give info for each profile included in the fitting
if true
    % k-path: W1-(20)-L21-(25)-G46-(29)-X75-(10)-U85-(1)-K86-(31)-G117
    profile(1).dir = '../profiles/GaN_zb_unstrained_0K_noSOC/'; %path to profile folder's location from this file's location
    profile(1).dft_file = 'dft.csv'; %filename containing DFT data
    profile(1).etb_data = ["GaN"]; %list of all .etb file needed for ETB calculation
    profile(1).dft_nkpoint = 117; %number of kpoints in DFT data
    profile(1).dft_vbo = 0.0; %value of DFT valence band edge
    profile(1).dft_band_list = [1:8]; %indices of bands in DFT data taken into fitting, in increasing order
    profile(1).etb_band_list = [1:8]; %indices of bands in calculated ETB to be used for fitting, must compatible with the `dft_band_list`
    profile(1).dft_top_vb = 4; %index of the top valence band in the dft.csv
    profile(1).etb_top_vb = 4; %index of the top valence band in the etb.csv
    %profile(1).coef2_list = [1:5]; %indices of bands in DFT data taken into fitting coef^2
    %profile(1).klist = [1:5:21, 26:5:46, 50:5:75, 80:5:85, 90:5:117]; %index of kpoints in DFT data taken into fitting
    profile(1).klist = [1:116];
    %profile(1).ene_window = [-5, 8]; %energy window for fitting (in eV), with respect to the DFT VBE
    profile(1).gauss = [10, 10; 1, 1]; % [amplitude_vb, amplitude_cb; width_vb, width_cb] for Gaussian band-weighting function (in eV)
    profile(1).bands_weight = [1, 10, 20, 20, 20, 2, 2, 2]; %weights for bands, length must be equal to that of `band_list`
    profile(1).coef2_weight = [1, 2, 3, 4, 5, 0, 0, 0]; %weights for coef^2, length must be equal to that of `band_list`
    profile(1).k_weight = [10, ones(1,19), 10, ones(1,24), 20, ones(1,28), 10, ones(1,9), 10, 1, ones(1,30)]; %weights for k-points
    profile(1).priority = 1;
    %profile(1).dft_bands = 'bands.csv'; %filename containing DFT bands
    %profile(1).vb_edge = 4; %in target file, the first `vb_edge` bands are valence: integer
    %profile(1).cb_edge = 5; %and from the `cb_edge` to the end are conduction: integer
    %profile(1).n_vband = 4; %number of valence bands taken into fitting (counted from VBE): integer
    %profile(1).n_cband = 6; %number of conduction bands taken into fitting (counted from CBE): integer
    %profile(1).n_cal_vband = 8; %number of valence bands actually calculated in TiberCAD simulation, not smaller than n_vband
    %profile(1).n_kpoint = 41; %number of k-points of each band in DFT target: integer
    %profile(1).crit_k_ind = [1, 21, 41]; %index of considered critical k-points to be paid more attention in `bands.dat.gnu` file: from `*.bands.in` file
    %profile(1).crit_k_name = ['L', 'G', 'X'];
    %profile(1).k_ind = [1, 11, 21, 26, 31, 36, 41]; %index of k-points in DFT target taken into fitting, e.g. all k-points or only critical ones
    %profile(1).focused_band = [3:10, 11];
end

%-----------------------------------
%OPTIMIZATION-CONTROLLING PARAMETERS
%-----------------------------------
%modifiable part
n_profile = 1; %number of profiles for fitting at the same time
n_objective = 2; %number of fitting objectives: 1 (only bands), or 2 (bands and coef2)
%tolerance = 0.0001; %tolerance for the cost value, in eV. Meaning: the averaged difference between ETB and DFT for each (E,k) point in bandtructure that you desire
%tol = 1e-4;             % Tolerance scale to compare two positions and keep the unique only
%algorithm = 'MOGA'; %choose one of the followings: 'SO' (Surrogate Optimization), 'PS' (Pattern Search), 'GA' (Genetic Algorithm), 'PSO' (Particle Swarm Optimizatipn), 'SA' (Simulated Annealing)
%hybrid = [];%'fgoalattain'; %use second optimizer after the first to refine the result. Only applied for 'PSO', 'GA', 'SA'. Built-in choices: [] (no use), 'fmincon' (recommended), 'patternsearch', 'fminunc', 'fminsearch'
parallel = true; %run in parallel or not (note: 'SA' does not support parallel mode)
%kpoint_priority = 10.0; %real number, = 1 to treat all kpoints equally, >1 to pay more attention to critical kpoints
%band_priority = 100.0; %real number, = 1 to treat all bands equally, >1 to pay more attention to focused bands (usually those near the gap)
%display = 'off'; %how matlab will print the progress to the terminal: 'iter', 'final' (default), 'off'
%guessed_point = transpose(importdata('selected+4.csv')); %bound for fitting range from the previous hydrostatic profile
weighting_method = 2; %0: uniform weight; 1: use band indices; 2: Gaussian of energy
bands_bound = 0.2; % cutoff for band cost
coef2_bound = 0.195; % cutoff for coef2 cost
mode = 0; % 0: start from scratch; 1: continue 

if mode == 0
    if isfile('../aux/estimated-point_stage1.csv')
        guessed_point = importdata('../aux/estimated-point_stage1.csv');
    else
        disp("ERROR: Cannot find file `../aux/estimated-point_stage1.csv`");
        return;
    end
elseif mode == 1
    if isfile('../aux/filtered_rep-point_stage2.csv')
        guessed_point = importdata('../aux/filtered_rep-point_stage2.csv');
    else
        disp("ERROR: Cannot find file `../aux/filtered_rep-point_stage2.csv`");
        return;
    end
else
    guessed_point = [];
end

%set options for the chose algorithm. See the Matlab documentation for the lists of these options.
%If not specififed, the options will be set to default.
options = optimoptions('gamultiobj', ...
    'PopulationSize', 3, ...  % Population size, 200
    'MaxGenerations', 3, ...  % Maximum number of generations, 200*n_etb_para, %'CrossoverFraction', 0.8, ...  % Fraction of population to replace in each generation, 0.8, %'EliteCount', 1, ...  % Number of elite individuals to keep for next generation, 5% of PopularSize
    'FunctionTolerance', 1e-5, ...  % Termination tolerance on the function value, 1e-4
    'OutputFcn', @outfunmoga, ... % Write out log after each iteration, %'PlotFcn', @gaplotbestf, ...  % Plot function to display the best fitness value in each generation
    'UseParallel', parallel, ... %'HybridFcn', hybrid, ...
    'Display', 'off', ...
    'MaxStallGenerations', 50, ... % 100, 
    'InitialPopulationMatrix', guessed_point, ...
    'ParetoFraction', 0.25); %Scalar from 0 through 1 specifying the fraction of individuals to keep on the first Pareto front, 0.35;

para_name = [...
    %"Delta_c_";... % representative for Delta_c + 4*Delta_correction_c_a
    %"Delta_a_";... % representative for Delta_a + 4*Delta_correction_a_c
    %"E_sc_";... %representative for E_sc + 4*I_sc*exp() + 4*O*exp() + offset --> keep fixed at 0.0
    "E_pc_";... % representative for E_pc + 4*I_pc*exp() + 4*O*exp() + offset
    "E_ec_";... % representative for E_ec + 4*I_ec*exp() + 4*O*exp() + offset
    "E_dc_";... % representative for E_dc + 4*I_dc*exp() + 4*O*exp() + offset
    "E_sa_";... % representative for E_sa + 4*I_sa*exp() + 4*O*exp() + offset
    "E_pa_";... % representative for E_pa + 4*I_pa*exp() + 4*O*exp() + offset
    "E_ea_";... % representative for E_ea + 4*I_ea*exp() + 4*O*exp() + offset
    %"E_da_";... % representative for E_da + 4*I_da*exp() + 4*O*exp() + offset
    "V_sss_";...
    "V_sps_";...
    "V_ses_";...
    "V_sds_";...
    "V_pss_";...
    "V_pps_";...
    "V_pes_";...
    "V_pds_";...
    "V_ess_";...
    "V_eps_";...
    "V_ees_";...
    "V_eds_";...
    %"V_dss_";...
    %"V_dps_";...
    %"V_des_";...
    %"V_dds_";...
    "V_ppp_";...
    "V_pdp_";...
    %"V_dpp_";...
    %"V_ddp_";...
    %"V_ddd_";...
    %"I_sc_a_";...
    %"I_pc_a_";...
    %"I_ec_a_";...
    %"I_dc_a_";...
    %"I_sa_c_";...
    %"I_pa_c_";...
    %"I_ea_c_";...
    %"I_da_c_";...
    %"delta_c_a_";...
    %"delta_a_c_",...
];

para_bound = [...
    %[   0.0246,   0.3000];... %Delta_c + 4*Delta_correction_c_a
    %[   0.0000,   0.0100];... %Delta_a + 4*Delta_correction_a_c
    %[  0.0000,   0.0000];... %E_sc + 4*I_sc*exp() + 4*O*exp() + offset = 0.0000
    [   5.0000,  15.0000];... %E_pc + 4*I_pc*exp() + 4*O*exp() + offset
    [  25.0000,  35.0000];... %E_ec + 4*I_ec*exp() + 4*O*exp() + offset
    [   7.0000,  27.0000];... %E_dc + 4*I_dc*exp() + 4*O*exp() + offset
    [ -16.0000,  -6.0000];... %E_sa + 4*I_sa*exp() + 4*O*exp() + offset
    [  -1.0000,   9.0000];... %E_pa + 4*I_pa*exp() + 4*O*exp() + offset
    [  15.0000,  35.0000];... %E_ea + 4*I_ea*exp() + 4*O*exp() + offset
    %[   7.0000,  27.0000];... %E_da + 4*I_da*exp() + 4*O*exp() + offset
    [ -10.0000,   0.0000];... %V_sss
    [   0.0000,  10.0000];... %V_sps
    [ -10.0000,   0.0000];... %V_ss*s
    [ -10.0000,   0.0000];... %V_sds
    [   0.0000,  10.0000];... %V_pss
    [   0.0000,  10.0000];... %V_pps
    [   0.0000,  10.0000];... %V_ps*s
    [ -10.0000,   0.0000];... %V_pds    
    [ -10.0000,   0.0000];... %V_s*ss
    [   0.0000,  10.0000];... %V_s*ps
    [ -10.0000,   0.0000];... %V_s*s*s
    [ -10.0000,   0.0000];... %V_s*ds
    %[ -10.0000,   0.0000];... %V_dss
    %[ -10.0000,   0.0000];... %V_dps
    %[ -10.0000,   0.0000];... %V_ds*s
    %[ -10.0000,   0.0000];... %V_dds
    [ -10.0000,   0.0000];... %V_ppp
    [   0.0000,  10.0000];... %V_pdp
    %[   0.0000,  10.0000];... %V_dpp
    %[   0.0000,  10.0000];... %V_ddp
    %[ -10.0000,   0.0000];... %V_ddd
    %[   0.0000,  20.0000];... %I_sc_a
    %[   0.0000,  20.0000];... %I_pc_a
    %[   0.0000,  20.0000];... %I_ec_a
    %[   0.0000,  20.0000];... %I_dc_a
    % 0.0000;... %I_sa_c
    % 0.0000;... %I_pa_c
    % 0.0000;... %I_ea_c
    % 0.0000;... %I_da_c
    % 0.0000;... %delta_lambda_ca
    % 0.0000;... %delta_lambda_ac
    ];

%=========
%MAIN PART
%=========
%---------
%CONSTANTS
%---------
if length(para_bound) ~= length(para_name)
    disp("ERROR: `para_bound` and `para_name` must have the same length.");
    return;
elseif size(guessed_point,1) > 0
    if length(para_bound) ~= size(guessed_point,2)
        disp("ERROR: `guessed_point` and `para_bound` must have the same length.");
        return;
    end
end

for p = 1:n_profile

    % self-check
    if length(profile(p).dft_band_list) ~= length(profile(p).etb_band_list)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `dft_bands_list` and `etb_band_list` must have the same length."));
        return;
    elseif length(profile(p).dft_band_list) ~= length(profile(p).bands_weight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `bands_weight` and `band_list` must have the same length."));
        return;
    elseif length(profile(p).dft_band_list) ~= length(profile(p).coef2_weight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `coef2_weight` and `band_list` must have the same length."));
        return;
    elseif length(profile(p).klist) ~= length(profile(p).k_weight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `k_weight` and `klist` must have the same length."));
        return;
    end
    %profile(p).n_band = profile(p).n_vband + profile(p).n_cband; %total number of considered bands
    n_k_fit = length(profile(p).klist); %number of kpoints taken into fitting
    n_band_fit = length(profile(p).etb_band_list); %number of bands taken into fitting

    %import dft bands
    dft_raw = importdata(strcat(profile(p).dir,profile(p).dft_file));
    idx = [];
    for b = profile(p).dft_band_list
        idx = [idx, (b-1)*profile(p).dft_nkpoint + profile(p).klist];
        if b == profile(p).dft_top_vb
            dft_vbe = max(dft_raw((b-1)*profile(p).dft_nkpoint + profile(p).klist, 2)); %assume VBE and CBE are given in dft.csv
            dft_cbe = min(dft_raw((b)*profile(p).dft_nkpoint + profile(p).klist, 2));
        end
    end
    profile(p).dft_bands = dft_raw(idx, 2); %value of VBE in dft.csv is assumed to be 0 eV
    profile(p).dft_coef2 = dft_raw(idx, 3:end);

    profile(p).weight.bands = ones(n_k_fit*n_band_fit, 1);
    profile(p).weight.coef2 = ones(n_k_fit*n_band_fit, 1);
    %for k = 1:length(p.crit_k_ind) %if take all kpoints into fitting
    %	weight(p.k_ind == p.crit_k_ind(k),:) = kpoint_priority;
    %end
    for k = 1:length(profile(p).klist)
        idx = [k:n_k_fit:(n_k_fit*(n_band_fit-1)+k)];
        profile(p).weight.bands(idx, 1) = profile(p).weight.bands(idx, 1) * profile(p).k_weight(k);
        profile(p).weight.coef2(idx, 1) = profile(p).weight.coef2(idx, 1) * profile(p).k_weight(k);
    end

    if weighting_method == 0
        w_bands = 1;
        w_coef2 = 1;
    elseif weighting_method == 1
        w_bands = [];
        w_coef2 = [];
        for band = 1:n_band_fit
            w_bands = [w_bands; ones(n_k_fit,1)*profile(p).bands_weight(band)];
            w_coef2 = [w_coef2; ones(n_k_fit,1)*profile(p).coef2_weight(band)];
        end
    elseif weighting_method == 2
        i_vbe = 0;
        for b = profile(p).dft_band_list
            if b <= profile(p).dft_top_vb
                i_vbe = i_vbe + n_k_fit;
            end
        end
        w_vb = profile(p).gauss(1,1) * normpdf(profile(p).dft_bands(1:i_vbe,1), dft_vbe, profile(p).gauss(2,1));
        w_cb = profile(p).gauss(1,2) * normpdf(profile(p).dft_bands(i_vbe+1:end,1), dft_cbe, profile(p).gauss(2,2));
        w_bands = [w_vb; w_cb];
        w_coef2 = w_bands;
    else
        disp("ERROR: `weighting_method` must be 0, 1, or 2.");
        return;
    end

    profile(p).weight.bands = profile(p).weight.bands .* w_bands;
    profile(p).weight.coef2 = profile(p).weight.coef2 .* w_coef2;
    %profile(p).weight.bands(:,1) = profile(p).weight.bands(:,1) .* profile(p).bands_weight(band);
    %profile(p).weight.coef2(:,1) = profile(p).weight.coef2(:,1) .* profile(p).coef2_weight(band);

    profile(p).weight.bands = profile(p).weight.bands/sum(profile(p).weight.bands, "all");
    profile(p).weight.coef2 = profile(p).weight.coef2/sum(profile(p).weight.coef2, "all");
    %profile(p).weight.bands = reshape(profile(p).weight.bands,[],1);
    %profile(p).weight.coef2 = reshape(profile(p).weight.coef2,[],1);
    %bottom = min(p.dft(1:p.n_k_ind)) - 0.5; %bottom value of the focused energy window
    %top = max(p.dft(p.n_vband*p.n_k_ind + 1:p.n_vband*p.n_k_ind + p.n_k_ind)) + 0.5; %top value of the focused energy window
    %weight(bottom < p.dft & p.dft < top) =  weight(bottom < p.dft & p.dft < top)*band_priority; %focus more on energies lying inside the energy window
    %weight = weight/sum(weight,"all");
    %for b=profile(p).band_list
    %    profile(p).dft_bands = [profile(p).dft_bands, dft_raw((b-1)*profile(p).dft_nkpoint + profile(p).klist, 2)];
    %end
    %profile(p).dft_bands = profile(p).dft_bands + profile(p).dft_vbo; %band structure is shifted to the value of DFT VBE

    %profile(p).dft_coef2 = dft_raw(profile(p).coef2_list)
    %for c = 1:size(dft_raw,2)-2 % iterate from column 3 to end in `dft.csv`: squared modules of s1, p1, d1, s2, p2, d2,...
    %    profile(p).dft_coef2{c} = [];
    %    for b = profile(p).coef2_list
    %        profile(p).dft_coef2{c} = [profile(p).dft_coef2{c}, dft_raw((b-1)*profile(p).dft_nkpoint + profile(p).klist, c+2)];
    %    end
    %end

    %import dft coef2
    %for k = 1:length(profile(p).crit_k_name)
    %    pdos = importdata(strcat(profile(p).dir, 'coef2_', profile(p).crit_k_name(k), '.csv'));
    %    pdos = [pdos; max(1.0 - sum(pdos,1), 0.0)]; %add a line of the sum of squared modulus for s* orbitals, put 0.0 if negative
    %    profile(p).dft_pdos{k} = transpose(pdos(:, [1:2, profile(p).vb_edge-profile(p).n_vband+3:profile(p).vb_edge+profile(p).n_cband]));
    %    %reduced coef2: sum of squares of same-type orbitals (sc, pc, dc, ec, sa, pa, da, ea)
    %    %to (sc+ec, all pc, all dc, sa+ea, all pa, all da)
    %    profile(p).dft_coef2{k}(:,1) = profile(p).dft_pdos{k}(:,1);
    %    profile(p).dft_coef2{k}(:,2) = sum(profile(p).dft_pdos{k}(:,2:4), 2);
    %    profile(p).dft_coef2{k}(:,3) = sum(profile(p).dft_pdos{k}(:,5:9), 2);
    %    profile(p).dft_coef2{k}(:,4) = profile(p).dft_pdos{k}(:,10);
    %    profile(p).dft_coef2{k}(:,5) = sum(profile(p).dft_pdos{k}(:,11:13), 2);
    %    profile(p).dft_coef2{k}(:,6) = sum(profile(p).dft_pdos{k}(:,14:18), 2);
    %    profile(p).dft_coef2{k}(:,7) = profile(p).dft_pdos{k}(:,19);
    %end

    %calculate the weights of the k-points in each profile according to focused energy window
    %profile(p).weight_bands = generate_weight(profile(p)); %for band cost
    %profile(p).weight_coef2 = profile(p).band_weight/sum(profile(p).band_weight,"all"); %for coef2 cost
    here = pwd;
    cd(profile(p).dir);
    system(strcat("rm -rf temp")); %remove all subfolders in the profile folder, if any
    system("mkdir temp");
    cd(here);
end
low_bound = transpose(para_bound(:,1)); %lower bounds
up_bound = transpose(para_bound(:,2)); %upper bounds
n_etb_para = size(para_bound,1); %number of ETB parameters: integer, first try: unstrained bulk

if ~isfile("../aux/rep-point_stage2.csv")
    system("touch ../aux/rep-point_stage2.csv");
end

if ~isfile("../aux/rep-costs_stage2.csv")
    system("touch ../aux/rep-costs_stage2.csv");
end

%---------------------------------------------
%DEFINE THE OBJECTIVE FUNCTION TO BE OPTIMIZED
%---------------------------------------------
object_fun = @(etb_set) cost_funs(n_profile, profile, para_name, etb_set, bands_bound, coef2_bound, n_objective);

%-------------------------------------
%OPTIMIZATION AND WRITE OUT THE RESULT
%-------------------------------------
if (parallel) %if choose parallel mode, then initiate it
    delete(gcp); %to ensure that no parallel pool exists
    parpool; %create a new parallel pool with a multicore processor, see: https://www.mathworks.com/help/gads/how-to-use-parallel-processing.html
end

%run algorithm
[opti_etb_set,cost_val,exitflag,~] = gamultiobj(object_fun, n_etb_para, [], [], [], [], low_bound, up_bound, [], options);

%export the results
%switch algorithm
%	case 'MOGA'
%		writematrix(transpose(opti_etb_set), 'pareto_points.csv');
%		writematrix(["#Band cost", "Coef2 cost"], 'pareto_costs.csv');
%		writematrix(cost_val, 'pareto_costs.csv', 'WriteMode','append');
%
%	otherwise
%		command_str = '';
%		for j=1:n_etb_para
%		    command_str = strcat(command_str,' -e "s/P',num2str(j),'_/',num2str(opti_etb_set(j),'%.6f'),'D0/g"');
%		end
%		total_cost_fun(n_profile, profile, n_etb_para, para_name, opti_etb_set, coef, true);
%		system(strcat('sed',command_str,' < sample.para > best.para'));
%end
time = toc;

%--------
%FINALIZE
%--------
fid = fopen('../aux/log_stage2.csv', 'a');
fprintf(fid, '%s %d\n', '#Exit flag: ', exitflag);
fprintf(fid, '%s %d\n', '#Min Bands Cost: ', min(cost_val(:,1)));
fprintf(fid, '%s %d\n', '#Min Coef2 Cost: ', min(cost_val(:,2)));
fprintf(fid, '%s %f\n', '#Run time (hours): ', time/3600);
fclose(fid);
disp('ALL FINISH!');
%exit %exit Matlab

%======================
%USER-DEFINED FUNCTIONS
%======================
%---------------------------------
function costs = cost_funs(n_profile, profile, para_name, etb_set, bands_bound, coef2_bound, n_objective)
    % Get current task information for unique ID
    %current_task = getCurrentTask();
    %if isempty(current_task)
    %    % If not in parallel mode, use random number
    %    temp_mate = strcat(num2str(randi(10)));
    %else
    %    % Use worker ID and process ID for uniqueness
    %    temp_mate = strcat(num2str(current_task.ID));
    %end
    
    % Ensure unique folder name
    %while (isfolder(strcat("../profiles/GaN_zb_unstrained_0K_noSOC/",temp_mate)))
    %    % If by chance the folder exists, append additional randomness
    %    temp_mate = strcat(temp_mate, '_', num2str(randi(10)));
    %end
    
    temp_mate = strcat(num2str(randi(1000000000))); %temporary material name for each instance of simulation, useful in case of parallel mode
    while (isfolder(strcat(profile(1).dir,temp_mate)))
        temp_mate = strcat(num2str(randi(1000000000)));
    end

    command_str = '';
    for j=1:length(etb_set)
        command_str = strcat(command_str,' -e "s/',para_name(j),'/',num2str(etb_set(j),'%.6f'),'D0/g"'); %create a command string for substituting multi-patterns
    end
    parent_dir = pwd;

    bands_cost = 0.0;
    coef2_cost = 0.0;
    priority = 0.0;
    for p = 1:n_profile
        cd(parent_dir);
        cd(strcat(profile(p).dir, "/temp/"));
        system(strcat("mkdir ", temp_mate));
        for f=1:length(profile(p).etb_data)
            system(strcat("sed",command_str," < ../../../datafiles/sample_",profile(p).etb_data(f),"_stage2.etb > ",temp_mate,"/",profile(p).etb_data(f),".etb")); %run sed to obtain the new temporary .etb file
        end
        full_klist = readmatrix("../kpoints_full");
        klist = full_klist(profile(p).klist,:);
        writematrix(klist,strcat(temp_mate,"/kpoints.csv"), 'Delimiter', ',');
        %run simulation
        %system(strcat('sed -e "s/@outdir@/out_',temp_mate,'/g" -e "s/@mate@/',temp_mate,'/g" < sample.tib > ',temp_mate,'.tib'));
        %system(strcat("cp sample.dat ",temp_mate,"/GaN.dat"));
        system(strcat("cp ../tb.upg ",temp_mate,"/tb.upg"));
        %system(strcat("cp kpoints_LkGkX ",temp_mate,"/kpoints_LkGkX"));
        system(strcat('sed -e "s|@temp@|',strcat(pwd,"/",temp_mate),'|g" < ../parameterizer.in > ',temp_mate,'/parameterizer.in'));
        %system(strcat("LD_LIBRARY_PATH= /home/alphan/TiberCAD/trunk_playground/bin/tibercad ",temp_mate,".tib")); %run TiberCAD
        %system(strcat("LD_LIBRARY_PATH= TIBERCADROOT=/home/alphan/TiberCAD/trunk_playground/ /home/alphan/TiberCAD/trunk_playground/x86_64-linux/bin/tibercad ",temp_mate,".tib")); %(Daniele: implementing Matthias' suggestion)
        cd(strcat(temp_mate));
        system(strcat("LD_LIBRARY_PATH= TIBERCADROOT=/home/alphan/TiberCAD/trunk_test/ /home/alphan/TiberCAD/trunk_test/src/core/modules/etb/libuptight/src/lib_uptight/PARAMETERIZER < parameterizer.in"));
        
        %import result and remove the temporary files and folder
        %if (isfile(strcat('out_',temp_mate,'/tb_dispersion.dat')))
        if (isfile('etb.csv'))
            %etb = importdata(strcat('out_',temp_mate,'/tb_dispersion.dat'));
            etb = importdata('etb.csv');
            idx = [];
            for b = profile(p).etb_band_list
                idx = [idx, (b-1)*length(klist) + [1:length(klist)]];
                if (b == profile(p).etb_top_vb)
                    idx_top_vb = (b-1)*length(klist) + [1:length(klist)];
                end
            end
            
            etb_bands = etb(idx, 2) - max(etb(idx_top_vb, 2));
            %selected_band_index = [profile(p).n_cal_vband+1:-1:profile(p).n_cal_vband-profile(p).n_vband+2, profile(p).n_cal_vband+2:profile(p).n_cal_vband+1+profile(p).n_cband];
            %etb = reshape(etb(:,2:end),1,[]) - max(etb(:,profile(p).n_cal_vband+1)); %band structure is shifted so that VBE = 0.0eV
            %etb = reshape(etb(:,2:end),1,[]);
            diff_bands = abs(etb_bands - profile(p).dft_bands).^2;
            %diff(isnan(diff)) = 0.0; %assign `diff` to 0 for points where DFT data is not available (i.e. not include these points into consideration)
            bands_cost = bands_cost + dot(profile(p).weight.bands, diff_bands)*profile(p).priority;

            if (n_objective == 2)
                etb_coef2 = etb(idx, 3:end);
                diff_coef2 = HellingerDist(etb_coef2, profile(p).dft_coef2);
                coef2_cost = coef2_cost + dot(profile(p).weight.coef2, diff_coef2)*profile(p).priority;
            end
            %for ao = 1:length(profile(p).dft_coef2)
            %    coef2 = importdata('pdos.dat');
            %    etb_coef2 = [];
            %    nbands = length(profile(p).band_list);
            %    nk = length(profile(p).klist);
            %    for e = 1:nbands
            %        etb_coef2 = [etb_coef2, coef2.data(e:nbands:nbands*nk, ao+2)];
            %    end
            %    diff_coef2 = abs(etb_coef2 - profile(p).dft_coef2{ao}); %implement Hellinger distance
            %    coef2_cost = coef2_cost + dot(profile(p).weight_coef2(:), diff_coef2(:).^2)*profile(p).priority;
            %end

            %for k = 1:length(profile(p).crit_k_name)
            %    etb_pdos = importdata(strcat('coef2_', profile(p).crit_k_name(k)));
            %    etb_coef2(:,1) = sum(etb_pdos(:,1:2), 2);
            %    etb_coef2(:,2) = sum(etb_pdos(:,3:5), 2);
            %    etb_coef2(:,3) = sum(etb_pdos(:,6:10), 2);
            %    etb_coef2(:,4) = sum(etb_pdos(:,11:12), 2);
            %    etb_coef2(:,5) = sum(etb_pdos(:,13:15), 2);
            %    etb_coef2(:,6) = sum(etb_pdos(:,16:20), 2);
            %    etb_coef2(:,7) = zeros(size(etb_pdos,1), 1);
            %    %etb_pdos_swap = etb_pdos; %swap order of two degeneracy states of opposite spins
            %    %for row = 1:2:size(etb_pdos,1)
            %    %    etb_pdos_swap(row,:) = etb_pdos(row+1,:);
            %    %    etb_pdos_swap(row+1,:) = etb_pdos(row,:);
            %    %end
            %    %diff_temp = sqrt(sum((etb_pdos - profile(p).dft_pdos{k}).^2, 2) / 2.0); %calculate Hellinger Distance for each eigenvector
            %    %diff_temp_swap = sqrt(sum((etb_pdos_swap - profile(p).dft_pdos{k}).^2, 2) / 2.0); %or for the swapped eigenvectors
            %    diff_temp = sqrt(sum((etb_coef2 - profile(p).dft_coef2{k}).^2, 2) / 2.0);
            %    %pdos_cost = pdos_cost +  min(dot(profile(p).weight_coef, diff_temp), dot(profile(p).weight_coef, diff_temp_swap))*profile(p).priority/numel(profile(p).crit_k_ind);
            %    pdos_cost = pdos_cost +  dot(profile(p).weight_coef2, diff_temp)*profile(p).priority/numel(profile(p).crit_k_ind);
            %end

        else
            disp(strcat("ERROR IN RUNNING TiberCAD FOR THE PROFILE ",num2str(p)));
            cd(parent_dir);
            return;
        end

	    %remove the temporary files of the previous run
        %system(strcat("rm -rf out_",temp_mate));
        %system(strcat("rm -r ",temp_mate,".*"));
        %system("rm -r *");
        %cd('..');
        %system(strcat("rm -rf ",temp_mate));

        priority = priority + profile(p).priority;
    end
    cd(parent_dir);
    costs = [bands_cost, coef2_cost]/priority;

    %bands_bound = 0.004;
    if (costs(1) > bands_bound) %to remove the unnecessary long tail of Pareto front in large band cost range
        costs(2) = costs(2)*10;
    end

    %coef2_bound = 0.07;
    if (costs(2) > coef2_bound) %to remove the unnecessary long tail of Pareto front in large coef cost range
        costs(1) = costs(1)*10;
    end

end

%----------------------------------------
% list of outfunc for different algorithm
%----------------------------------------
function [state,options,optchanged] = outfunmoga(options,state,flag)
    switch flag
        case 'init'
            % Initialize the CSV file
            fid = fopen('../aux/log_stage2.csv', 'w');
            fprintf(fid, '#Multi-Objective Genetic Algorithm\n#Generation, Min Bands Cost, Min Coef2 Cost\n');
            fclose(fid);
        case 'iter'
        	min_band = min(state.Score(:,1));
        	min_pdos = min(state.Score(:,2));

            isDominated = IsDominated(state.Score);
            rep_point = importdata('../aux/rep-point_stage2.csv');
            rep_point = [rep_point; state.Population(~isDominated,:)];
            rep_costs = importdata('../aux/rep-costs_stage2.csv');
            rep_costs = [rep_costs; state.Score(~isDominated,:)];

            isDominated = IsDominated(rep_costs);
            rep_point = rep_point(~isDominated,:);
            rep_costs = rep_costs(~isDominated,:);

            [rep_point, i_old, ~] = uniquetol(rep_point, 1e-4, 'ByRows', true); %keep only unique positions
            rep_costs = rep_costs(i_old,:); %new rep_costs after remove duplicated elements

            writematrix(rep_point, '../aux/rep-point_stage2.csv');
            writematrix(rep_costs, '../aux/rep-costs_stage2.csv');
            %writematrix(state.Score, 'costs.csv');
            %writematrix(state.Population, 'population.csv');

            fid = fopen('../aux/log_stage2.csv', 'a');
            fprintf(fid, '%d, %f, %f\n', state.Generation, min_band, min_pdos);
            fclose(fid);
        case 'done'
            fid = fopen('../aux/log_stage2.csv', 'a');
            fprintf(fid, '\n');
            fclose(fid);
    end
    optchanged = false;
end

%------------------
function isDominated = IsDominated(costs)

    nPop = size(costs,1);
    isDominated = false(nPop,1);
    
    for i = 1:nPop-1
        for j = i+1:nPop
            
            if (all(costs(i,:) <= costs(j,:)) && any(costs(i,:) < costs(j,:))) %if j is dominated by i
                isDominated(j) = true;
            end
            
            if (all(costs(j,:) <= costs(i,:)) && any(costs(j,:) < costs(i,:))) %if i is dominated by j
                isDominated(i) = true;
            end
            
        end
    end

end

%------------------
function dist = HellingerDist(x, y)
    dist = [];
    for i = 1:size(x, 1)
        %dist = [dist; sqrt(sum((sqrt(x(i,:)) - sqrt(y(i,:))).^2, 2) / 2.0)];
        dist = [dist; sqrt(sum((x(i,:) - y(i,:)).^2, 2) / 2.0)];
    end
end
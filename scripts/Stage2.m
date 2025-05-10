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
    profile(1).dir = '../profiles/GaN_zb_unstrained_300K_noSOC/'; %path to profile folder's location from this file's location
    profile(1).dft_file = 'dft.csv'; %filename containing DFT data
    profile(1).etb_data = ["GaN"]; %list of all .etb file needed for ETB calculation
    profile(1).dft_nkpoint = 117; %number of kpoints in DFT data
    profile(1).dft_vbo = 0.0; %value of DFT valence band edge
    profile(1).dft_band_list = [1:6]; %indices of bands in DFT data taken into fitting, in increasing order
    profile(1).etb_band_list = [1:6]; %indices of bands in calculated ETB to be used for fitting, must compatible with the `dft_band_list`
    profile(1).dft_top_vb = 4; %index of the top valence band in the dft.csv
    profile(1).etb_top_vb = 4; %index of the top valence band in the etb.csv
    profile(1).klist = [1:5:46, 50:5:85, 89:6:113];
    profile(1).gauss = [6, 4]; % [gaussian-sigma_vb, sigma_cb] for Gaussian band-weighting function (in eV)
    profile(1).bands_bweight = [5, 3, 3, 3, 3, 4]; %factor weights for bands_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(1).coef2_bweight = [0, 1, 1, 1, 0, 0]; %factor weights for coef2_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(1).bands_kweight = [15, 3, 20, 3, 30, 7, 2, 2, 5, 50, 5, 8, 15, 7, 5, 12, 5, 15, 8, 10, 10, 4, 3]; %weights for bands_cost in terms of k-index, must have the same length as `klist`
    profile(1).coef2_kweight = [15, 5, 5, 8, 20, 7, 7, 15, 15, 35, 15, 10, 5, 3, 5, 20, 10, 20, 5, 5, 8, 7, 5]; %weights for coef2_cost in terms of k-index, must have the same length as `klist`
    profile(1).priority = 1;
    %============================================

    % k-path: W1-(20)-L21-(25)-G46-(29)-X75-(10)-U85-(1)-K86-(31)-G117
    profile(2).dir = '../profiles/InN_zb_unstrained_300K_noSOC/'; %path to profile folder's location from this file's location
    profile(2).dft_file = 'dft.csv'; %filename containing DFT data
    profile(2).etb_data = ["InN"]; %list of all .etb file needed for ETB calculation
    profile(2).dft_nkpoint = 117; %number of kpoints in DFT data
    profile(2).dft_vbo = 0.0; %value of DFT valence band edge
    profile(2).dft_band_list = [1:6]; %indices of bands in DFT data taken into fitting, in increasing order
    profile(2).etb_band_list = [1:6]; %indices of bands in calculated ETB to be used for fitting, must compatible with the `dft_band_list`
    profile(2).dft_top_vb = 4; %index of the top valence band in the dft.csv
    profile(2).etb_top_vb = 4; %index of the top valence band in the etb.csv
    profile(2).klist = [1, 4, 8, 15, 21:5:46, 52:6:70, 75:5:85, 92:6:110];
    profile(2).gauss = [5.5, 6]; % [gaussian-sigma_vb, sigma_cb] for Gaussian band-weighting function (in eV)
    profile(2).bands_bweight = [5, 3, 3, 3, 3, 5]; %factor weights for bands_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(2).coef2_bweight = [0, 1, 1.5, 1.5, 1, 0]; %factor weights for coef2_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(2).bands_kweight = [30, 10, 15, 5, 25, 10, 10, 15, 20, 50, 20, 15, 15, 5, 25, 5, 20, 20, 20, 8, 4]; %weights for bands_cost in terms of k-index, must have the same length as `klist`
    profile(2).coef2_kweight = [15, 3, 0, 15, 25, 15, 15, 15, 15, 50, 15, 10, 15, 5, 20, 5, 20, 0, 10, 10, 10]; %weights for coef2_cost in terms of k-index, must have the same length as `klist`
    profile(2).priority = 1;

    % k-path: W1-(20)-L21-(25)-G46-(29)-X75-(10)-U85-(1)-K86-(31)-G117
    profile(3).dir = '../profiles/AlN_zb_unstrained_300K_noSOC/'; %path to profile folder's location from this file's location
    profile(3).dft_file = 'dft.csv'; %filename containing DFT data
    profile(3).etb_data = ["AlN"]; %list of all .etb file needed for ETB calculation
    profile(3).dft_nkpoint = 117; %number of kpoints in DFT data
    profile(3).dft_vbo = 0.0; %value of DFT valence band edge
    profile(3).dft_band_list = [1:6]; %indices of bands in DFT data taken into fitting, in increasing order
    profile(3).etb_band_list = [1:6]; %indices of bands in calculated ETB to be used for fitting, must compatible with the `dft_band_list`
    profile(3).dft_top_vb = 4; %index of the top valence band in the dft.csv
    profile(3).etb_top_vb = 4; %index of the top valence band in the etb.csv
    profile(3).klist = [1, 5:5:15, 21:5:46, 51:4:75, 80, 85, 88, 92, 97:7:111];
    profile(3).gauss = [5.0, 4.4]; % [gaussian-sigma_vb, sigma_cb] for Gaussian band-weighting function (in eV)
    profile(3).bands_bweight = [5, 3, 3, 3, 3, 4]; %factor weights for bands_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(3).coef2_bweight = [0, 1, 1, 1, 1, 0]; %factor weights for coef2_cost in terms of band index, length must be equal to that of `dft_band_list`
    profile(3).bands_kweight = [25, 5, 1, 15, 25, 20, 20, 20, 20, 100, 20, 10, 20, 20, 15, 10, 20, 15, 10, 50, 20, 20, 10, 3]; %weights for bands_cost in terms of k-index, must have the same length as `klist`
    profile(3).coef2_kweight = [10, 5, 5, 0, 10, 5, 10, 20, 30, 50, 30, 20, 10, 8, 8, 8, 10, 5, 15, 2, 5, 2, 5, 5]; %weights for coef2_cost in terms of k-index, must have the same length as `klist`
    profile(3).priority = 1;
end%

%-----------------------------------
%OPTIMIZATION-CONTROLLING PARAMETERS
%-----------------------------------
%modifiable part
profile_list = [2]; %list of profiles to be used for fitting
n_objective = 2; %number of fitting objectives: 1 (only bands), or 2 (bands and coef2)
parallel = true; %run in parallel or not (note: 'SA' does not support parallel mode)
bands_bound = 0.034; % cutoff for band cost
coef2_bound = 0.062; % cutoff for coef2 cost
mode = 1; % 0: no initial guess; 1: continue from a previous cycle; else: use estimation from Stage 1

if mode == 0
    guessed_point = [];
elseif mode == 1
    if isfile('../aux/filtered_rep-point_stage2.csv')
        guessed_point = importdata('../aux/filtered_rep-point_stage2.csv');
    else
        disp("ERROR: Cannot find file `../aux/filtered_rep-point_stage2.csv`");
        return;
    end
else
    if isfile('../aux/estimated.csv')
        guessed_point = importdata('../aux/estimated.csv');
    else
        disp("ERROR: Cannot find file `../aux/estimated.csv`");
        return;
    end
end

%set options for the chose algorithm. See the Matlab documentation for the lists of these options.
%If not specififed, the options will be set to default.
options = optimoptions('gamultiobj', ...
    'PopulationSize', 123, ...  % Population size, 200
    'MaxGenerations', Inf, ...  % Maximum number of generations, 200*n_etb_para, %'CrossoverFraction', 0.8, ...  % Fraction of population to replace in each generation, 0.8, %'EliteCount', 1, ...  % Number of elite individuals to keep for next generation, 5% of PopularSize
    'FunctionTolerance', 1e-5, ...  % Termination tolerance on the function value, 1e-4
    'OutputFcn', @outfunmoga, ... % Write out log after each iteration, %'PlotFcn', @gaplotbestf, ...  % Plot function to display the best fitness value in each generation
    'UseParallel', parallel, ... %'HybridFcn', hybrid, ...
    'Display', 'off', ...
    'MaxStallGenerations', 40, ... % 100, 
    'InitialPopulationMatrix', guessed_point, ...
    'ParetoFraction', 0.4); %Scalar from 0 through 1 specifying the fraction of individuals to keep on the first Pareto front, 0.35;

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
    "E_da_";... % representative for E_da + 4*I_da*exp() + 4*O*exp() + offset
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
    "V_dss_";...
    "V_dps_";...
    "V_des_";...
    "V_dds_";...
    "V_ppp_";...
    "V_pdp_";...
    "V_dpp_";...
    "V_ddp_";...
    "V_ddd_";...
];

para_bound = [...
    %[   0.0246,   0.3000];... %Delta_c + 4*Delta_correction_c_a
    %[   0.0000,   0.0100];... %Delta_a + 4*Delta_correction_a_c
    %[  0.0000,   0.0000];... %E_sc + 4*I_sc*exp() + 4*O*exp() + offset = 0.0000
    [   8.0000,  13.0000];... %E_pc + 4*I_pc*exp() + 4*O*exp() + offset
    [  18.0000,  23.0000];... %E_ec + 4*I_ec*exp() + 4*O*exp() + offset
    [  22.0000,  27.0000];... %E_dc + 4*I_dc*exp() + 4*O*exp() + offset
    [ -14.0000, -10.0000];... %E_sa + 4*I_sa*exp() + 4*O*exp() + offset
    [  -5.0000,   0.0000];... %E_pa + 4*I_pa*exp() + 4*O*exp() + offset
    [  21.0000,  25.0000];... %E_ea + 4*I_ea*exp() + 4*O*exp() + offset
    [  14.0000,  18.0000];... %E_da + 4*I_da*exp() + 4*O*exp() + offset
    [ -20.0000,   0.0000];... %V_sss
    [   0.0000,  20.0000];... %V_sps
    [ -20.0000,   0.0000];... %V_ss*s
    [ -20.0000,   0.0000];... %V_sds
    [   0.0000,  20.0000];... %V_pss
    [   0.0000,  20.0000];... %V_pps
    [   0.0000,  20.0000];... %V_ps*s
    [ -20.0000,   0.0000];... %V_pds    
    [ -20.0000,   0.0000];... %V_s*ss
    [   0.0000,  20.0000];... %V_s*ps
    [ -20.0000,   0.0000];... %V_s*s*s
    [ -20.0000,   0.0000];... %V_s*ds
    [ -20.0000,   0.0000];... %V_dss
    [ -20.0000,   0.0000];... %V_dps
    [ -20.0000,   0.0000];... %V_ds*s
    [ -20.0000,   0.0000];... %V_dds
    [ -20.0000,   0.0000];... %V_ppp
    [   0.0000,  20.0000];... %V_pdp
    [   0.0000,  20.0000];... %V_dpp
    [   0.0000,  20.0000];... %V_ddp
    [ -20.0000,   0.0000];... %V_ddd
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

for p = profile_list

    % self-check
    if length(profile(p).dft_band_list) ~= length(profile(p).etb_band_list)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `dft_bands_list` and `etb_band_list` must have the same length."));
        return;
    elseif length(profile(p).dft_band_list) ~= length(profile(p).bands_bweight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `bands_weight` and `band_list` must have the same length."));
        return;
    elseif length(profile(p).dft_band_list) ~= length(profile(p).coef2_bweight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `coef2_weight` and `band_list` must have the same length."));
        return;
    elseif length(profile(p).klist) ~= length(profile(p).bands_kweight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `bands_kweight` and `klist` must have the same length."));
        return;
    elseif length(profile(p).klist) ~= length(profile(p).coef2_kweight)
        disp(strcat("ERROR IN THE PROFILE ",num2str(p),": `coef2_kweight` and `klist` must have the same length."));
        return;
    end
    
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

    % calculate the weights in terms of k-points
    profile(p).weight.bands = ones(n_k_fit*n_band_fit, 1);
    profile(p).weight.coef2 = ones(n_k_fit*n_band_fit, 1);
    
    for k = 1:length(profile(p).klist)
        idx = [k:n_k_fit:(n_k_fit*(n_band_fit-1)+k)];
        profile(p).weight.bands(idx, 1) = profile(p).weight.bands(idx, 1) * profile(p).bands_kweight(k);
        profile(p).weight.coef2(idx, 1) = profile(p).weight.coef2(idx, 1) * profile(p).coef2_kweight(k);
    end

    % calculate the weights in terms of energy
    i_vbe = 0;
    for b = profile(p).dft_band_list
        if b <= profile(p).dft_top_vb
            i_vbe = i_vbe + n_k_fit;
        end
    end
    w_vb = normpdf(profile(p).dft_bands(1:i_vbe,1), dft_vbe, profile(p).gauss(1,1));
    w_cb = normpdf(profile(p).dft_bands(i_vbe+1:end,1), dft_cbe, profile(p).gauss(1,2));
    w_bands = [w_vb; w_cb];
    w_coef2 = w_bands;

    fac_bands = [];
    fac_coef2 = [];
    for band = 1:n_band_fit
        fac_bands = [fac_bands; ones(n_k_fit,1)*profile(p).bands_bweight(band)];
        fac_coef2 = [fac_coef2; ones(n_k_fit,1)*profile(p).coef2_bweight(band)];
    end

    w_bands = w_bands .* fac_bands;
    w_coef2 = w_coef2 .* fac_coef2;

    % calculate the net weights
    profile(p).weight.bands = profile(p).weight.bands .* w_bands;
    profile(p).weight.coef2 = profile(p).weight.coef2 .* w_coef2;

    profile(p).weight.bands = profile(p).weight.bands/sum(profile(p).weight.bands, "all");
    profile(p).weight.coef2 = profile(p).weight.coef2/sum(profile(p).weight.coef2, "all");

    here = pwd;
    cd(profile(p).dir);
    if ~isfolder("temp")
        system("mkdir temp");
    end
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
object_fun = @(etb_set) cost_funs(profile_list, profile, para_name, etb_set, bands_bound, coef2_bound, n_objective);

%-------------------------------------
%OPTIMIZATION AND WRITE OUT THE RESULT
%-------------------------------------
if (parallel) %if choose parallel mode, then initiate it
    delete(gcp); %to ensure that no parallel pool exists
    parpool; %create a new parallel pool with a multicore processor, see: https://www.mathworks.com/help/gads/how-to-use-parallel-processing.html
end

%run algorithm
[opti_etb_set,cost_val,exitflag,~] = gamultiobj(object_fun, n_etb_para, [], [], [], [], low_bound, up_bound, [], options);

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
function costs = cost_funs(profile_list, profile, para_name, etb_set, bands_bound, coef2_bound, n_objective)
    
    % Ensure unique folder name
    %while (isfolder(strcat("../profiles/GaN_zb_unstrained_0K_noSOC/",temp_mate)))
    %    % If by chance the folder exists, append additional randomness
    %    temp_mate = strcat(temp_mate, '_', num2str(randi(10)));
    %end
    
    temp_mate = strcat(num2str(randi(1000000000))); %temporary material name for each instance of simulation, useful in case of parallel mode
    trial = 0;
    while ((isfolder(strcat(profile(1).dir,temp_mate))) && trial < 1000)
        temp_mate = strcat(num2str(randi(1000000000)));
        trial = trial + 1;
    end

    command_str = '';
    for j=1:length(etb_set)
        command_str = strcat(command_str,' -e "s/',para_name(j),'/',num2str(etb_set(j),'%.6f'),'D0/g"'); %create a command string for substituting multi-patterns
    end
    parent_dir = pwd;

    bands_cost = 0.0;
    coef2_cost = 0.0;
    priority = 0.0;
    for p = profile_list
        cd(parent_dir);
        cd(strcat(profile(p).dir, "/temp/"));
        system(strcat("mkdir ", temp_mate));
        for f=1:length(profile(p).etb_data)
            system(strcat("sed",command_str," < ../../../datafiles/sample_",profile(p).etb_data(f),"_stage2.etb > ",temp_mate,"/",profile(p).etb_data(f),".etb")); %run sed to obtain the new temporary .etb file
        end
        full_klist = readmatrix("../kpoints");
        klist = full_klist(profile(p).klist,:);
        writematrix(klist,strcat(temp_mate,"/kpoints.csv"), 'Delimiter', ',');

        %run simulation
        system(strcat("cp ../tb.upg ",temp_mate,"/tb.upg"));
        system(strcat('sed -e "s|@temp@|',strcat(pwd,"/",temp_mate),'|g" < ../parameterizer.in > ',temp_mate,'/parameterizer.in'));
        cd(strcat(temp_mate));
        system(strcat("LD_LIBRARY_PATH= TIBERCADROOT=/home/alphan/TiberCAD/trunk_test/ /home/alphan/TiberCAD/trunk_test/src/core/modules/etb/libuptight/src/lib_uptight/PARAMETERIZER < parameterizer.in"));
        
        %import result and remove the temporary files and folder
        if (isfile('etb.csv'))
            etb = importdata('etb.csv');
            idx = [];
            for b = profile(p).etb_band_list
                idx = [idx, (b-1)*length(klist) + [1:length(klist)]];
                if (b == profile(p).etb_top_vb)
                    idx_top_vb = (b-1)*length(klist) + [1:length(klist)];
                end
            end
            
            etb_bands = etb(idx, 2) - max(etb(idx_top_vb, 2));
            diff_bands = abs(etb_bands - profile(p).dft_bands).^2;
            bands_cost = bands_cost + dot(profile(p).weight.bands, diff_bands)*profile(p).priority;

            if (n_objective == 2)
                etb_coef2 = etb(idx, 3:end);
                diff_coef2 = HellingerDist(etb_coef2, profile(p).dft_coef2);
                coef2_cost = coef2_cost + dot(profile(p).weight.coef2, diff_coef2)*profile(p).priority;
            end

        else
            disp(strcat("ERROR IN RUNNING TiberCAD FOR THE PROFILE ",num2str(p)));
            cd(parent_dir);
            return;
        end

        cd('..');
        system(strcat("rm -rf ",temp_mate));

        priority = priority + profile(p).priority;
    end
    cd(parent_dir);
    costs = [bands_cost, coef2_cost]/priority;

    % Redefine the costs to remove the unnecessary long tail of Pareto front in large band cost range
    if (costs(1) > bands_bound)
        costs(2) = costs(2)*10;
    end

    if (costs(2) > coef2_bound)
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
        %dist = [dist; sqrt(sum((sqrt(x(i,:)) - sqrt(y(i,:))).^2, 2) / 2.0)]; % the Hellinger formula by definition
        dist = [dist; sqrt(sum((x(i,:) - y(i,:)).^2, 2) / 2.0)]; % modified formula to tolerate the "negative orbital characters"
    end
end

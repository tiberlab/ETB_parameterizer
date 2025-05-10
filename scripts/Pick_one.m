function Pick_one(profile_list, material, stage, series, save, date, index)

    %INPUTS: profile_list: list of profile indices, e.g. [1,2]
    %        material: name of the material, e.g. "GaN"
    %        stage: stage number in string, e.g. "2"
    %        series: series number in string, e.g. "1"
    %        save: save number in string, e.g. "07"
    %        date: today in string, e.g. "26Apr2025"
    %        index: index of the chosen point in the repository in the order of ascending band cost

    %cost = importdata('pareto_costs.csv');
    %cost = cost.data;
    cost = importdata(strcat('../archives/',material,'/Stage',stage,'/series',series,'_save',save,'_',date,'/rep-costs_stage2.csv'));
    [sorted_band_cost, i_particle] = sort(cost(:,1));

    disp( strcat("BAND COST: ", num2str(sorted_band_cost(index)) ) );
    disp( strcat("COEF COST: ", num2str(cost(i_particle(index),2)) ) );
    %if (sorted_band_cost(index) > band_threshold)
    %	disp('BAND COST IS OVER THRESHOLD!');
    %	return;
    %end
    
    profile(1).dir = '../profiles/GaN_zb_unstrained_300K_noSOC/';
    profile(1).etb_data = ["GaN"];

    profile(2).dir = '../profiles/InN_zb_unstrained_300K_noSOC/';
    profile(2).etb_data = ["InN"];

    profile(3).dir = '../profiles/AlN_zb_unstrained_300K_noSOC/';
    profile(3).etb_data = ["AlN"];
    
    %particle = importdata('pareto_points.csv');
    particle = transpose(importdata(strcat('../archives/',material,'/Stage',stage,'/series',series,'_save',save,'_',date,'/rep-point_stage2.csv')));
    para = particle(:,i_particle(index));
    n_etb_para = length(para);
    %writematrix(para, 'selected0.csv');

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

    
    command_str = '';
    for j=1:n_etb_para
        command_str = strcat(command_str,' -e "s/',para_name(j),'/',num2str(para(j),'%.6f'),'D0/g"');
    end
    
    parent_dir = pwd;
    for p = profile_list
        cd(parent_dir);
        cd(profile(p).dir);
        %system("rm -rf ./*/");
        system(strcat("mkdir select2"));
        for f=1:length(profile(p).etb_data)
            system(strcat("sed",command_str," < ../../datafiles/sample_",profile(p).etb_data(f),"_stage2.etb > select2/",profile(p).etb_data(f),".etb")); %run sed to obtain the new temporary .etb file
        end

        %run simulation
        %system(strcat('sed -e "s/@outdir@/out_',temp_mate,'/g" -e "s/@mate@/',temp_mate,'/g" < sample.tib > ',temp_mate,'.tib'));
        %system(strcat("cp sample.dat ",temp_mate,"/GaN.dat"));
        system(strcat("cp tb.upg select2/tb.upg"));
        system(strcat("cp kpoints select2/kpoints.csv"));
        system(strcat('sed -e "s|@temp@|',strcat(pwd,"/select2"),'|g" < parameterizer.in > select2/parameterizer.in'));
        %system(strcat("LD_LIBRARY_PATH= /home/alphan/TiberCAD/trunk_playground/bin/tibercad ",temp_mate,".tib")); %run TiberCAD
        %system(strcat("LD_LIBRARY_PATH= TIBERCADROOT=/home/alphan/TiberCAD/trunk_playground/ /home/alphan/TiberCAD/trunk_playground/x86_64-linux/bin/tibercad ",temp_mate,".tib")); %(Daniele: implementing Matthias' suggestion)
        cd("select2/");
        system(strcat("LD_LIBRARY_PATH= TIBERCADROOT=/home/alphan/TiberCAD/trunk_test/ /home/alphan/TiberCAD/trunk_test/src/core/modules/etb/libuptight/src/lib_uptight/PARAMETERIZER < parameterizer.in"));
        
        %import result and remove the temporary files and folder
        %if (isfile(strcat('out_',temp_mate,'/tb_dispersion.dat')))
        if (~isfile('etb.csv'))
            disp(strcat("ERROR IN RUNNING TiberCAD FOR THE PROFILE ",num2str(p)));
            cd(parent_dir);
            return;
        end
    end

    cd(parent_dir);
    disp( strcat("BAND COST: ", num2str(sorted_band_cost(index)) ) );
    disp( strcat("COEF COST: ", num2str(cost(i_particle(index),2)) ) );

    return;
end
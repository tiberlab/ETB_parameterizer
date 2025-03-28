rep_point = importdata('../temp/rep-point_stage2.csv');
rep_costs = importdata('../temp/rep-costs_stage2.csv');

band_filter = 0.19;
coef_filter = 0.188;
step = 6;

filtered_ind = all([rep_costs(:,1) < band_filter, rep_costs(:,2) < coef_filter], 2);

filtered_point = rep_point(filtered_ind,:);
filtered_costs = rep_costs(filtered_ind,:);

filtered_point = filtered_point(1:step:end,:);
filtered_costs = filtered_costs(1:step:end,:)

writematrix(filtered_point, '../temp/filtered_rep-point_stage2.csv');

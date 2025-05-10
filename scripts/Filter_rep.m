function Filter_rep(band_filter, coef_filter, step)
% Filter_rep filters representative points based on band and coefficient costs.
%
%   Inputs:
%       band_filter: Maximum allowed band cost.
%       coef_filter: Maximum allowed coefficient cost.
%       step: Step size for downsampling the filtered points.

    % Load data
    rep_point = importdata('../aux/rep-point_stage2.csv');
    rep_costs = importdata('../aux/rep-costs_stage2.csv');

    % Filter points based on cost thresholds
    filtered_ind = all([rep_costs(:,1) < band_filter, rep_costs(:,2) < coef_filter], 2);

    filtered_point = rep_point(filtered_ind,:);
    filtered_costs = rep_costs(filtered_ind,:);

    % Downsample the filtered points
    filtered_point = filtered_point(1:step:end,:);

    % Display the number of filtered points
    num_filtered_points = size(filtered_point,1);
    disp(['Number of filtered points: ', num2str(num_filtered_points)])

    % Write the filtered points to a file
    writematrix(filtered_point, '../aux/filtered_rep-point_stage2.csv');

end

%% PD fits

addpath(genpath('D:\causality\Causality Perception'));

% load PD and Control data
pd = table2array(readtable("patient.xlsx", "ReadVariableNames", true, 'sheet', 2));
ct = table2array(readtable("controls.xlsx", "ReadVariableNames", true, 'sheet', 1));
% define the overlaps

x = pd(:,1);
x = x*-1;

a_pd = zeros(15, 1);
b_pd = zeros(15, 1);
hf_pd = zeros(15, 1);
amp_pd = zeros(15, 1);

for ii = 1:15;
    temp_dat = [x,pd(:,ii+1)];
    [a_pd(ii), b_pd(ii), hf_pd(ii), amp_pd(ii)] = fitting_log(temp_dat);
    xlim([-1 4]);
    ylim([0 1.05]);
    hold on
end

x = pd(:,1); 
x = x*-1
ys = pd(:,2:end);   

long_x = repmat(x, size(ys, 2), 1);
long_y = reshape(ys, [], 1);

dat = [long_x, long_y];

[~,~,~,~,plt] = fitting_log(dat)
set(plt, 'LineWidth', 2, 'Color', 'r');

%% Control fits
a_ct = zeros(15, 1);
b_ct = zeros(15, 1);
hf_ct = zeros(15, 1);
amp_ct = zeros(15, 1);

figure;

for ii = 1:15
    temp_dat = [x,ct(:,ii+1)];
    [a_ct(ii), b_ct(ii), hf_ct(ii), amp_ct(ii)] = fitting_log(temp_dat);
    xlim([-1 4]);
    ylim([0 1.05]);
end

x = ct(:,1);
x = x*-1
ys = ct(:,2:end);


long_x = repmat(x, size(ys, 2), 1);
long_y = reshape(ys, [], 1);

dat = [long_x, long_y];

[~,~,~,~,plt] = fitting_log(dat)
set(plt, 'LineWidth', 2, 'Color', 'b'); 

%% Stat test between pd and ct
[h, p, ci, stat] = ttest2(a_ct,a_pd,"Vartype","unequal")

cohend = computeCohen_d(a_ct,a_pd,"paired")

%% Error Bars
x = pd(:,1); 
x = x*-1
ys = pd(:,2:end);    

long_x = repmat(x, size(ys, 2), 1);
long_y = reshape(ys, [], 1);

unique_x1 = unique(long_x);

y_mean1 = zeros(size(unique_x1));
y_sem1 = zeros(size(unique_x1));

% Calculate mean and SEM for each unique x
for i = 1:length(unique_x1)
    y_values = long_y(long_x == unique_x(i));
    y_mean1(i) = mean(y_values);
    y_sem1(i) = std(y_values) / sqrt(length(y_values));
end

% Plot scatter with error bars
figure;
errorbar(unique_x1, y_mean1, y_sem1, 'o', 'LineWidth', 1.2, 'MarkerSize', 5, 'CapSize', 5,'Color','r');
ylim([0 1.05])
xlim([-0.5 4])
hold on

% controls
x = ct(:,1); 
x = x*-1
ys = ct(:,2:end); 

long_x = repmat(x, size(ys, 2), 1);
long_y = reshape(ys, [], 1);

unique_x1 = unique(long_x);

y_mean1 = zeros(size(unique_x1));
y_sem1 = zeros(size(unique_x1));

% Calculate mean and SEM for each unique x
for i = 1:length(unique_x1)
    y_values = long_y(long_x == unique_x(i));
    y_mean1(i) = mean(y_values);
    y_sem1(i) = std(y_values) / sqrt(length(y_values));
end

% Plot scatter with error bars
errorbar(unique_x1, y_mean1, y_sem1, 'o', 'LineWidth', 1.2, 'MarkerSize', 5, 'CapSize', 5,'Color','b');
ylim([0 1.05])
xlim([-0.5 4])
hold off

%% Box plots

% Combine data
as  = [a_pd,  a_ct];  
hfs = [hf_pd, hf_ct];  

n = size(as,1);

% Jitter for scatter points
jitter = -0.2 + (0.4)*rand(n, 2); 

x_l = ones(n,2);
x_l(:,2) = 2;
x_l2 = x_l + jitter;

% Plotting
figure;
subplot(1,2,1)
boxchart(x_l(:), as(:), 'BoxFaceColor','w','BoxEdgeColor','k','BoxWidth',0.4)
hold on
scatter(x_l2(:,1), as(:,1), 30, 'r', 'filled','MarkerEdgeColor','k')
scatter(x_l2(:,2), as(:,2), 30, 'b', 'filled','MarkerEdgeColor','k')
ylabel('Slope')
xticks([1 2])
xticklabels({'PD','Control'})
hold off

subplot(1,2,2)
boxchart(x_l(:), hfs(:), 'BoxFaceColor','w','BoxEdgeColor','k','BoxWidth',0.4)
hold on
scatter(x_l2(:,1), hfs(:,1), 30, 'r', 'filled')
scatter(x_l2(:,2), hfs(:,2), 30, 'b', 'filled')
ylabel('PSE')
xticks([1 2])
xticklabels({'PD','Control'})
hold off


%% reaction times
pd_rt = table2array(readtable("patient.xlsx", "ReadVariableNames", true, 'sheet', 3));
pd_rt=pd_rt(1:110,:)

ct_rt = table2array(readtable("controls.xlsx", "ReadVariableNames", true, 'sheet', 2));
ct_rt=ct_rt(1:110,:)


% Extract unique x-values and subjects data
x_values = unique(pd_rt(:, 1)); 
num_x = length(x_values);       
num_subjects = size(pd_rt, 2) - 1; 

% Initialize means and SEM for each x-value
pd_means = zeros(num_x, 1);
pd_sem = zeros(num_x, 1);
ct_means = zeros(num_x, 1);
ct_sem = zeros(num_x, 1);

% Loop through unique x-values and calculate mean and SEM
for i = 1:num_x
    pd_idx = pd_rt(:, 1) == x_values(i);
    ct_idx = ct_rt(:, 1) == x_values(i);

    pd_data = pd_rt(pd_idx, 2:end);
    ct_data = ct_rt(ct_idx, 2:end);

    pd_means(i) = mean(pd_data(:));
    pd_sem(i) = std(pd_data(:)) / sqrt(num_subjects);

    ct_means(i) = mean(ct_data(:));
    ct_sem(i) = std(ct_data(:)) / sqrt(num_subjects);
end

xx = [1 2 3 4 5 6 7 8 9 10 11]

% Create the bar plot
figure;
hold on;
bar(xx - 0.2, pd_means, 0.3, 'FaceColor', 'r', 'EdgeColor', 'none');
bar(xx + 0.2, ct_means, 0.3, 'FaceColor', 'b', 'EdgeColor', 'none');
errorbar(xx - 0.2, pd_means, pd_sem, 'k.', 'LineWidth', 0.5);
errorbar(xx + 0.2, ct_means, ct_sem, 'k.', 'LineWidth', 0.5);
xlabel('Overlaps');
ylabel('Mean Response Reaction Time');
legend({'PD', 'Control'}, 'Location', 'Best');
hold off;

%% ----------------rt analysis-------------------------------
% Combine data
overlaps_pd = pd_rt(:, 1); 
overlaps_ct = ct_rt(:, 1); 

responses_pd = pd_rt(:, 2:end); 
responses_ct = ct_rt(:, 2:end); 

% Get number of subjects
nSubj_pd = size(responses_pd, 2);
nSubj_ct = size(responses_ct, 2);

% Create group labels
group_pd = repmat({'PD'}, numel(responses_pd), 1);
group_ct = repmat({'Control'}, numel(responses_ct), 1);

% Create subject labels
subjects_pd = repelem(1:nSubj_pd, size(responses_pd, 1))';
subjects_ct = repelem(1:nSubj_ct, size(responses_ct, 1))';

% Create long format data
overlap_long = [repmat(overlaps_pd, nSubj_pd, 1); repmat(overlaps_ct, nSubj_ct, 1)];
response_long = [responses_pd(:); responses_ct(:)];
group_long = [group_pd; group_ct];
subjects_long = [subjects_pd; subjects_ct];

% Create table
dataTable = table(overlap_long, response_long, group_long, subjects_long, ...
    'VariableNames', {'Overlap', 'Response', 'Group', 'Subject'});

% ✅ Convert to categorical
dataTable.Overlap = categorical(dataTable.Overlap);  
dataTable.Group = categorical(dataTable.Group);
dataTable.Subject = categorical(dataTable.Subject);

% Fit the LME model
lme_model = fitlme(dataTable, ...
    'Response ~ Group + Overlap + (1 | Subject)', ...
    'FitMethod', 'REML');

disp(lme_model);

anovaTable = anova(lme_model, 'DFMethod', 'Residual');
disp(anovaTable);

VarRandom = sum(var(lme_model.randomEffects)); % Variance from random effects
VarResidual = lme_model.MSE;                   % Residual variance (error)

FittedValues = fitted(lme_model, 'Conditional', false);
VarFixed = var(FittedValues);                  % Variance explained by fixed effects

% Total variance
VarTotal = VarFixed + VarRandom + VarResidual;

% Marginal R-squared (fixed effects only)
R2_Marginal = VarFixed / VarTotal;

% Conditional R-squared (fixed + random effects)
R2_Conditional = (VarFixed + VarRandom) / VarTotal;

fprintf('Marginal R^2: %.4f\n', R2_Marginal);
fprintf('Conditional R^2: %.4f\n', R2_Conditional);

%% ----------------------repeated mesure ANOVA----------------------------

% === PD group ===

unique_overlaps_pd = unique(pd_rt(:,1),'stable');
nOverlap = numel(unique_overlaps_pd);

nSubj_pd = size(pd_rt,2)-1;
nReps_pd = sum(pd_rt(:,1)==unique_overlaps_pd(1));

pd_data_matrix = zeros(nSubj_pd, nOverlap);

% Loop over subjects
for s = 1:nSubj_pd
    subj_data = pd_rt(:,s+1);
    for o = 1:nOverlap
        this_overlap = unique_overlaps_pd(o);
        idx = pd_rt(:,1) == this_overlap;
        pd_data_matrix(s,o) = mean(subj_data(idx));
    end
end

clean_overlap_labels = strrep(string(unique_overlaps_pd'), '.', '_');
varnames = strcat('Overlap_', clean_overlap_labels);

pd_table = array2table(pd_data_matrix, 'VariableNames', varnames);
pd_table.Subject = strcat('PD_', string(1:nSubj_pd))';
pd_table.Group = repmat({'PD'}, nSubj_pd, 1);

% === Control group ===

unique_overlaps_ct = unique(ct_rt(:,1),'stable');

if ~isequal(unique_overlaps_pd, unique_overlaps_ct)
    warning('Overlap levels differ between PD and Control groups.');
end

nSubj_ct = size(ct_rt,2)-1;
nReps_ct = sum(ct_rt(:,1)==unique_overlaps_ct(1));

ct_data_matrix = zeros(nSubj_ct, nOverlap);

for s = 1:nSubj_ct
    subj_data = ct_rt(:,s+1);
    
    for o = 1:nOverlap
        this_overlap = unique_overlaps_ct(o);
        idx = ct_rt(:,1) == this_overlap;
        
        ct_data_matrix(s,o) = mean(subj_data(idx));
    end
end

ct_table = array2table(ct_data_matrix, 'VariableNames', varnames);
ct_table.Subject = strcat('CT_', string(1:nSubj_ct))';
ct_table.Group = repmat({'Control'}, nSubj_ct, 1);

anova_table = [pd_table; ct_table];

anova_table.Subject = categorical(anova_table.Subject);
anova_table.Group = categorical(anova_table.Group);

WithinDesign = table(categorical(unique_overlaps_pd), 'VariableNames', {'Overlap'});
WithinDesign.Properties.RowNames = varnames;

% === Fit repeated measures ANOVA model ===

rm = fitrm(anova_table, ...
    strcat(varnames{1}, '-', varnames{end}, ' ~ Group'), ...
    'WithinDesign', WithinDesign, 'WithinModel', 'Overlap');

ranova_results = ranova(rm, 'WithinModel', 'Overlap');
disp(ranova_results);


%% overlaps
overpd = table2array(readtable("C:\Users\hadic\Downloads\overlaps\overlapsp.xlsx", "ReadVariableNames", true, 'sheet', 2));
overct = table2array(readtable("C:\Users\hadic\Downloads\overlaps\overlapc.xlsx", "ReadVariableNames", true, 'sheet', 2));

x = overpd(:,1)'
y1 = overpd(:,2:end)
aa = 5

y2 = overct(:,2:end)
mean_y1 = mean(y1,2)'
sem_y1 = std(y1')/sqrt(length(y1(1,:)))
mean_y2 = mean(y2,2)'
sem_y2 = std(y2')/sqrt(length(y2(1,:)))

figure
% Construct X and Y coordinates for the patch
x_patch = [x, fliplr(x)]; 
y1_patch = [mean_y1 + sem_y1, fliplr(mean_y1 - sem_y1)]; 
y2_patch = [mean_y2 + sem_y2, fliplr(mean_y2 - sem_y2)];
% Plot the SEM as a properly shaped patch (shaded area)
patch(x_patch, y1_patch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
patch(x_patch, y2_patch, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(x, mean_y1, 'r', 'LineWidth', 1.2); 
plot(x, mean_y2, 'b', 'LineWidth', 1.2);
ylim([0 4])
xlabel('Overlaps'); 
ylabel('Normalized Overlap Perception Score');
hold off;


figure
hold on;
errorbar(x, mean_y1, sem_y1, 'ro', 'LineWidth', 2, 'CapSize', 20); 
errorbar(x, mean_y2, sem_y2, 'bo', 'LineWidth', 2, 'CapSize', 20); 
ylim([0 3]);
xlim([-10 110])
xlabel('Overlaps'); 
ylabel('Normalized Overlap Perception Score');
hold off;

%% R-square on overlap

% Combine data
overlaps_pd = overpd(:, 1); 
overlaps_ct = overct(:, 1); 

responses_pd = overpd(:, 2:end); 
responses_ct = overct(:, 2:end);

nSubj_pd = size(responses_pd, 2);
nSubj_ct = size(responses_ct, 2);

group_pd = repmat({'PD'}, numel(responses_pd), 1);
group_ct = repmat({'Control'}, numel(responses_ct), 1);

subjects_pd = repelem(1:nSubj_pd, size(responses_pd, 1))';
subjects_ct = repelem(1:nSubj_ct, size(responses_ct, 1))';

overlap_long = [repmat(overlaps_pd, nSubj_pd, 1); repmat(overlaps_ct, nSubj_ct, 1)];
response_long = [responses_pd(:); responses_ct(:)];
group_long = [group_pd; group_ct];
subjects_long = [subjects_pd; subjects_ct];

dataTable = table(overlap_long, response_long, group_long, subjects_long, ...
    'VariableNames', {'Overlap', 'Response', 'Group', 'Subject'});

dataTable.Group = categorical(dataTable.Group);
dataTable.Subject = categorical(dataTable.Subject);

% Fit the LME model
lme_model = fitlme(dataTable, ...
    'Response ~ Group + Overlap + (1|Subject)', ...
    'FitMethod', 'REML');

disp(lme_model);

% Extract variance components
VarRandom = sum(var(lme_model.randomEffects)); % Variance from random effects
VarResidual = lme_model.MSE; % Residual variance (error)

% Compute fitted values for fixed effects
FittedValues = fitted(lme_model, 'Conditional', false);
VarFixed = var(FittedValues); % Variance explained by fixed effects

% Total variance
VarTotal = VarFixed + VarRandom + VarResidual;

% Marginal R-squared (fixed effects only)
R2_Marginal = VarFixed / VarTotal;

% Conditional R-squared (fixed + random effects)
R2_Conditional = (VarFixed + VarRandom) / VarTotal;

fprintf('Marginal R^2: %.4f\n', R2_Marginal);
fprintf('Conditional R^2: %.4f\n', R2_Conditional);

%% -----------------repeated measure on overlap----------------------------

% === PD group ===

% Extract unique overlaps
unique_overlaps_pd = overpd(:,1);
nOverlap = numel(unique_overlaps_pd);

nSubj_pd = size(overpd,2)-1;

pd_data_matrix = overpd(:,2:end)'; % rows = subjects, cols = overlaps

clean_overlap_labels = strrep(string(unique_overlaps_pd'), '.', '_');
varnames = strcat('Overlap_', clean_overlap_labels);

pd_table = array2table(pd_data_matrix, 'VariableNames', varnames);
pd_table.Subject = strcat('PD_', string(1:nSubj_pd))';
pd_table.Group = repmat({'PD'}, nSubj_pd, 1);

% === Control group ===

unique_overlaps_ct = overct(:,1);

if ~isequal(unique_overlaps_pd, unique_overlaps_ct)
    warning('Overlap levels differ between PD and Control groups.');
end

nSubj_ct = size(overct,2)-1;

ct_data_matrix = overct(:,2:end)';

ct_table = array2table(ct_data_matrix, 'VariableNames', varnames);
ct_table.Subject = strcat('CT_', string(1:nSubj_ct))';
ct_table.Group = repmat({'Control'}, nSubj_ct, 1);

% === Combine both groups ===

anova_table = [pd_table; ct_table];

anova_table.Subject = categorical(anova_table.Subject);
anova_table.Group = categorical(anova_table.Group);

% === Define WithinDesign table ===

overlap_levels = categorical(unique_overlaps_pd);
WithinDesign = table(overlap_levels, 'VariableNames', {'Overlap'});
WithinDesign.Properties.RowNames = varnames;

% === Fit repeated measures ANOVA model ===

rm = fitrm(anova_table, ...
    strcat(varnames{1}, '-', varnames{end}, ' ~ Group'), ...
    'WithinDesign', WithinDesign, 'WithinModel', 'Overlap');

ranova_results = ranova(rm, 'WithinModel', 'Overlap');
disp(ranova_results);

%% speeds

data = readtable("D:\causality\ratios_final3.csv"); 
data = data(data.speed ~= 0.5, :);

% Define unique levels
unique_speeds = unique(data.speed);
unique_overlaps = unique(data.overlap);

% Set consistent colors for overlaps
cmap = lines(length(unique_overlaps));

% Prepare mean and SEM matrices
mean_matrix = zeros(length(unique_speeds), length(unique_overlaps));
sem_matrix = zeros(length(unique_speeds), length(unique_overlaps));

for s = 1:length(unique_speeds)
    speed_level = unique_speeds(s);
    for o = 1:length(unique_overlaps)
        overlap_level = unique_overlaps(o);
        
        subset = data(data.speed == speed_level & ...
                      data.overlap == overlap_level, :);
        
        mean_matrix(s, o) = mean(subset.ratio);
        sem_matrix(s, o) = std(subset.ratio) / sqrt(height(subset));
    end
end

% Create figure
figure;
bh = bar(mean_matrix, 'grouped');
hold on;

for o = 1:length(bh)
    bh(o).FaceColor = cmap(o, :);
end

[numGroups, numBars] = size(mean_matrix);
groupWidth = min(0.8, numBars/(numBars + 1.5));
for o = 1:numBars
    x = (1:numGroups) - groupWidth/2 + (2*o-1) * groupWidth / (2*numBars);
    errorbar(x, mean_matrix(:, o), sem_matrix(:, o), 'k.', 'LineWidth', 1.2);
end

xticks(1:length(unique_speeds));
xticklabels(string(unique_speeds));
xlabel('Speed Level');
ylabel('Causality Ratio');
ylim([0 1]);
title('All Subjects');
legend(arrayfun(@(x) sprintf('Overlap %.1f', x), unique_overlaps, 'UniformOutput', false), ...
       'Location', 'northwest');
hold off;

%% Mixed ANOVA between speeds
data = readtable("D:\causality\ratios_final3.csv"); 
data = data(data.speed ~= 0.5, :);

% Ensure categorical variables
data.subject = categorical(data.subject);
data.speed = categorical(data.speed);
data.overlap = categorical(data.overlap);

lme = fitlme(data, ...
    'ratio ~ speed * overlap + (1|subject)');

fprintf('\n--- LME Results ---\n');
disp(lme)

fprintf('\n--- Mixed ANOVA ---\n');
disp(anova(lme)); 

% Extract variance components
var_fixed = var(predict(lme));
var_random = lme.covarianceParameters{1,1}; % random intercept variance
var_resid = lme.MSE;                        % residual variance

% Compute pseudo-R²
r2_marginal = var_fixed / (var_fixed + var_random + var_resid);
r2_conditional = (var_fixed + var_random) / (var_fixed + var_random + var_resid);

fprintf('Goodness of Fit (R²):\n');
fprintf('  Marginal R² (fixed effects):     %.4f\n', r2_marginal);
fprintf('  Conditional R² (fixed + random): %.4f\n', r2_conditional);
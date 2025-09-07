function plot_models_evaluation()
% Simplified version: Compare multiple models' ROC curves and metrics

% ===== 1. Model files configuration =====
files = {
    'score_matrix', 'model1';
    'score_matrix', 'model2'; ...

};

% ===== 2. Initialize =====
mean_fpr = linspace(0, 1, 100);
nModels = size(files, 1);
results = cell(nModels, 6); % Model, AUC, AUC_STD, PR_AUC, F1, Accuracy

% ===== 3. Plot setup =====
figure('Position', [100, 100, 780, 620]); 
hold on;
colors = lines(nModels);
legends = cell(nModels, 1);

% ===== 4. Process each model =====
for i = 1:nModels
    filename = files{i,1};
    label = files{i,2};
    
    % Load data
    data = csvread(filename);
    y_true = logical(data(:,1));
    scores_matrix = data(:, 2:end);
    
    % Calculate metrics for each round
    rounds = size(scores_matrix, 2);
    aucs = zeros(rounds, 1);
    pr_aucs = zeros(rounds, 1);
    f1s = zeros(rounds, 1);
    accs = zeros(rounds, 1);
    tprs_all = zeros(rounds, length(mean_fpr));
    
    for j = 1:rounds
        scores = scores_matrix(:, j);
        
        % ROC curve
        [fpr, tpr, ~, auc_val] = perfcurve(y_true, scores, 1);
        aucs(j) = auc_val;
        
        % PR curve
        [~, ~, ~, pr_auc_val] = perfcurve(y_true, scores, 1, 'xCrit', 'reca', 'yCrit', 'prec');
        pr_aucs(j) = pr_auc_val;
        
        % Best F1 and Accuracy (simplified threshold search)
        thresholds = linspace(0, 1, 101);
        f1_scores = zeros(size(thresholds));
        acc_scores = zeros(size(thresholds));
        
        for k = 1:length(thresholds)
            y_pred = scores >= thresholds(k);
            tp = sum(y_pred & y_true);
            fp = sum(y_pred & ~y_true);
            fn = sum(~y_pred & y_true);
            tn = sum(~y_pred & ~y_true);
            
            precision = tp / (tp + fp + eps);
            recall = tp / (tp + fn + eps);
            f1_scores(k) = 2 * precision * recall / (precision + recall + eps);
            acc_scores(k) = (tp + tn) / length(y_true);
        end
        
        f1s(j) = max(f1_scores);
        accs(j) = max(acc_scores);
        
        % Interpolate TPR for average ROC
        [fpr_unique, ia] = unique(fpr);
        tpr_unique = tpr(ia);
        tprs_all(j, :) = interp1(fpr_unique, tpr_unique, mean_fpr, 'linear', 'extrap');
    end
    
    % Calculate averages
    mean_auc = mean(aucs);
    std_auc = std(aucs);
    mean_pr_auc = mean(pr_aucs);
    mean_f1 = mean(f1s);
    mean_acc = mean(accs);
    mean_tpr = mean(tprs_all, 1);
    
    % Plot average ROC curve
    plot(mean_fpr, mean_tpr, 'LineWidth', 2, 'Color', colors(i,:));
    legends{i} = sprintf('%s (AUC = %.3f±%.3f)', label, mean_auc, std_auc);
    
    % Store results
    results(i, :) = {label, round(mean_auc,4), round(std_auc,4), ...
                     round(mean_pr_auc,4), round(mean_f1,4), round(mean_acc,4)};
    
    % Print results
    fprintf('%s: AUC=%.3f±%.3f, PR-AUC=%.3f, F1=%.3f, Acc=%.3f\n', ...
            label, mean_auc, std_auc, mean_pr_auc, mean_f1, mean_acc);
end

% ===== 5. Finalize plot =====
plot([0 1], [0 1], 'k--'); % Diagonal line
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('Average ROC Curves Comparison');
legend(legends, 'Location', 'SouthEast');
grid on;

% ===== 6. Save results =====
results_table = cell2table(results, 'VariableNames', ...
    {'Model','AUC','AUC_STD','PR_AUC','F1','Accuracy'});
writetable(results_table, 'model_evaluation_results.xlsx');
end
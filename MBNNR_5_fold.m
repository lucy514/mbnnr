% ======= Main Program Entry =======
A = load('data\Protein_Disease_adj.txt'); % Load protein-disease association matrix
nd = max(A(:,1));
nv = max(A(:,2));
[pp, qq] = size(A); % pp = 905

% Load interaction matrix and fixed negative sample file
interaction_matrix = xlsread('data\Protein_Disease_Associations.xlsx');
[row, col] = size(interaction_matrix);
number = 100; % Define number of rounds

% Initialize result storage matrix

all_results = zeros(1810, 101);
all_results(:,1) = [ones(905,1); zeros(905,1)]; % First column is label

% BNNR parameters
maxiter = 300;
alpha = 1;
beta = 1;
tol1 = 2e-3;
tol2 = 1e-5;

array1 = xlsread('data\positive_samples_100rounds.xlsx');

% Load negative sample file
fixed_negative_samples = xlsread('data\negative_samples_100rounds.xlsx');

tic;
for k = 1:100
    fprintf('Round %d\n', k);
    disp(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    array = array1(k,:);
    
    % Get all positive sample pairs
    [pos_row, pos_col] = find(interaction_matrix == 1);
    pos_pairs = [pos_row, pos_col];
    
    % Get all negative sample pairs
    [neg_row, neg_col] = find(interaction_matrix == 0);
    neg_pairs = [neg_row, neg_col];
    
    % Use negative samples (read negative sample indices for round k from file)
    fixed_neg_indices = fixed_negative_samples(k, :);
    selected_neg_pairs = [neg_row(fixed_neg_indices), neg_col(fixed_neg_indices)];
    
    % Initialize prediction scores for current round
    current_round_scores = zeros(1810, 1);
    
    % Five-fold cross validation
    for circle = 0:4
        nu = 181;
        interaction_matrix1 = interaction_matrix;
        
        % Get test indices for current fold
        if circle < 4
            new_array = array(1, 1 + circle * nu : (circle + 1) * nu);
        else
            new_array = array(1, 1 + circle * nu : end);
        end
        
        % Update training matrix
        for j = 1:length(new_array)
            o = A(new_array(j), 1);
            l = A(new_array(j), 2);
            interaction_matrix1(o, l) = 0;
        end
        
        % Dynamic Gaussian similarity calculation
        protein_gauss = similarity_protein(interaction_matrix1);
        disease_gauss = similarity_disease(interaction_matrix1);
        
        % Integration function to combine similarity matrices
        [protein_integration_similarity,disease_integration_similarity] = integration_protein_disease_similarity(protein_gauss, disease_gauss);
        
        
        % BNNR with masking
        Wdd = protein_integration_similarity;
        Wvv = disease_integration_similarity;
        Wvd = interaction_matrix1';
        [dn, dr] = size(Wvd);
        T = [Wdd, Wvd'; Wvd, Wvv];
        [t1, t2] = size(T);
        trIndex = double(T ~= 0);  
        [WW, ~] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
        M_recovery = (WW((end - dn + 1):end, 1:dr))';
        
        A1 = M_recovery;
        A_Opposite = (1 - interaction_matrix);
        A11 = A1 .* A_Opposite;
        col_sum = sum(A11);
        A11_norm = A11 ./ col_sum;
        A_M = interaction_matrix .* A1;
        predict_score_matrix = A11_norm + A_M;
        
        % Save prediction scores for current fold
        for j = 1:length(new_array)
            q = A(new_array(j), 1);
            w = A(new_array(j), 2);
            score = predict_score_matrix(q, w);
            current_round_scores(j + circle * nu) = score;
        end
        
        % Save negative sample prediction scores in the last fold
        if circle == 4
            for i = 1:size(selected_neg_pairs, 1)
                q = selected_neg_pairs(i, 1);
                w = selected_neg_pairs(i, 2);
                score = predict_score_matrix(q, w);
                current_round_scores(905 + i) = score; % Start storing negative sample scores from position 906
            end
        end
    end
    
    % Store current round results in total result matrix
    all_results(:, k+1) = current_round_scores;
end

toc;

% Create table headers
headers = cell(1, 101);
headers{1} = 'Label';
for i = 1:100
    headers{i+1} = ['Round' num2str(i)];
end

% Convert results to table
T = array2table(all_results, 'VariableNames', headers);

% Save results
writetable(T, 'prediction_scores_bnnr_hybrid_1-100.xlsx');



% Load protein number to Uniprot ID mapping table
[~, ~, raw] = xlsread('data\Protein_Numbers.xlsx');
protein_ids_all = cell2mat(raw(2:end,1));
uniprot_ids_all = raw(2:end,2);           

disease_id = 2; % Fill in the disease number you are interested in here
interaction_matrix = xlsread('data\Protein_Disease_Associations.xlsx');
protein_gauss = similarity_protein(interaction_matrix);
disease_gauss = similarity_disease(interaction_matrix);
        
% Integration function to combine similarity matrices
[protein_integration_similarity,disease_integration_similarity] = integration_protein_disease_similarity(protein_gauss, disease_gauss);
Wdd = protein_integration_similarity;
Wvv = disease_integration_similarity;
Wvd = interaction_matrix';
T = [Wdd, Wvd'; Wvd, Wvv];
trIndex = double(T ~= 0);

alpha = 1; 
beta = 1; 
tol1 = 2e-3; 
tol2 = 1e-5;
 maxiter = 300;
[WW, ~] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);

[dn, dr] = size(Wvd);
M_recovery = (WW((end - dn + 1):end, 1:dr))';
A1 = M_recovery;
A_Opposite = (1 - interaction_matrix);
A11 = A1 .* A_Opposite;
col_sum = sum(A11);
A11_norm = A11 ./ col_sum;
A_M = interaction_matrix .* A1;
predict_score_matrix = A11_norm + A_M;

% Find all negative samples (unassociated proteins) for this disease
neg_proteins = find(interaction_matrix(:, disease_id) == 0);

% Extract scores for these proteins
neg_scores = [];
for i = 1:length(neg_proteins)
    protein_id = neg_proteins(i);
    score = predict_score_matrix(protein_id, disease_id);
    neg_scores = [neg_scores; protein_id, score];
end

% Sort by score in descending order and take the top 30
neg_scores_sorted = sortrows(neg_scores, -2);
top30 = neg_scores_sorted(1:min(30, size(neg_scores_sorted,1)), :);

% Find UniprotID corresponding to protein numbers
protein_ids = top30(:,1);
top30_uniprot_ids = cell(size(protein_ids));
for i = 1:length(protein_ids)
    idx = find(protein_ids_all == protein_ids(i));
    if ~isempty(idx)
        top30_uniprot_ids{i} = uniprot_ids_all{idx};
    else
        top30_uniprot_ids{i} = '';
    end
end

% Output
Top30_table = table(protein_ids, top30_uniprot_ids, top30(:,2), 'VariableNames', {'ProteinID', 'UniprotID', 'Score'});
writetable(Top30_table, sprintf('Disease%d_Top30_negative_proteins_with_uniprot.xlsx', disease_id));
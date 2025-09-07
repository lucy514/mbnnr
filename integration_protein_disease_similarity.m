function [ protein_similarity,disease_similarity ] = integration_protein_disease_similarity(protein_gauss, disease_gauss)
%INTEGRATION_PROTEIN_DISEASE_SIMILARITY Integrate pre-computed and dynamic similarity matrices
%   [protein_similarity, disease_similarity] = INTEGRATION_PROTEIN_DISEASE_SIMILARITY(protein_gauss, disease_gauss)
%   integrates pre-computed sequence/semantic similarity with dynamic Gaussian similarity
%
%   Inputs:
%       protein_gauss - dynamic Gaussian similarity matrix for proteins
%       disease_gauss - dynamic Gaussian similarity matrix for diseases
%
%   Outputs:
%       protein_similarity - integrated protein similarity matrix
%       disease_similarity - integrated disease similarity matrix

% Note: Previous versions used weight matrices for weighted combination
% Current version: simple addition of similarity matrices

% Load pre-computed similarity matrices
ss = xlsread('data\Protein_Sequence_Similarity_Matrix.xlsx'); % Protein sequence similarity matrix
vs = xlsread('data\Disease_Semantic_Similarity_Matrix.xlsx'); % Disease semantic similarity matrix

% Integrate similarities by simple addition
protein_similarity = ss + protein_gauss;     % Protein integrated similarity matrix
disease_similarity = vs + disease_gauss;     % Disease integrated similarity matrix

end

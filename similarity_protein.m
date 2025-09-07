function KD=similarity_protein(interaction)
%SIMILARITY_PROTEIN Calculate Gaussian similarity matrix for proteins
%   KD = SIMILARITY_PROTEIN(interaction) calculates the Gaussian kernel
%   similarity matrix between proteins based on their interaction patterns
%
%   Inputs:
%       interaction - interaction matrix where rows represent proteins
%
%   Outputs:
%       KD - protein-protein similarity matrix

[nd,~]=size(interaction);
for i=1:nd
    sh(i)=norm(interaction(i,:))^2;  % calculate gamma for Gaussian kernel calculation
end
gamad=nd/sum(sh');
for i=1:nd
    for j=1:nd
        kh(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);  % calculate Gaussian kernel similarity
    end
end
KD=kh;

% Note: This function calculates protein-protein similarity using Gaussian kernel
% based on interaction patterns to avoid data leakage in cross-validation
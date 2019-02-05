function [dG0,reverse_ddG0] = legendretransformR(S,kegg_pKa,pH,T,I,dG0_prime)

R = 8.31e-3; % kJ/mol/K

fprintf('Performing reverse transform for reaction Gibbs energies\n');
% Calculate the reverse transform for all reactions  
% From Noor et al. (2013)
reverse_ddG0 = zeros(size(S, 2), 1);
I(isnan(I)) = 0.25; % default ionic strength is 0.25M
%pMg(isnan(pMg)) = 14; % default pMg is 14
for i = 1:size(S, 2) % for each reaction in S
    inds = find(S(:, i));
    reaction_ddG0s = zeros(length(inds), 1);
    for j = 1:length(inds)
        diss = kegg_pKa(inds(j));
        dG0s = cumsum(-[0, diag(diss.pKas, 1)'] * R * T(i) * log(10));
        dG0s = dG0s - dG0s(diss.majorMSpH7);
        pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)];
        reaction_ddG0s(j) = Transform(pseudoisomers, pH(i), I(i), T(i));
    end
    reverse_ddG0(i) = S(inds, i)' * reaction_ddG0s;
end

dG0 = dG0_prime - reverse_ddG0;

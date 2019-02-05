function [dG0,reverse_ddG0s] = legendretransformF(m,kegg_pKa,cids_m,pH,T,I,dG0_prime)

R = 8.31e-3; % kJ/mol/K
fprintf('Performing reverse transform for formation energies\n');
% Calculate the reverse transform for all formation energies 
I(isnan(I)) = 0.25; % default ionic strength is 0.25M
reverse_ddG0s = zeros(m, 1);
for j = 1:(m)
         % find the diss table from the pKa data
         k = find(strcmp(cids_m(j),{kegg_pKa.cid}));
            if ~isempty(k)
                diss = kegg_pKa(k);
            end
        
        if isempty(diss)
            continue;
        end
        
        dG0s = cumsum(-[0, diag(diss.pKas, 1)'] * R * T(j) * log(10));
        dG0s = dG0s - dG0s(diss.majorMSpH7);
        pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)];
        reverse_ddG0s(j) = Transform(pseudoisomers, pH(j), I(j), T(j));
end

dG0 = dG0_prime - reverse_ddG0s;

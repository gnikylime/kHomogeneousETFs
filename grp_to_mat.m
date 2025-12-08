function mat_ind = grp_to_mat(group,grp_ind)
% grp_to_mat Converting elements in fin abel grp to Matlab indices

% INPUTS:
% group: vector of natural numbers [n_1, ..., n_k] representing a finite 
% % abelian group G = Z_(n_1) x ... x Z_(n_k)
% grp_ind: vector of k integers representing an element of G, with the ith
% % entry in the range 0, ..., n_i - 1.

% OUTPUT:
% mat_ind: integer in the range 1, ..., n_1*n_2*...*n_k corresponding to number
% % of group element in lexicographic ordering

% Author: Emily J King
% https://www.math.colostate.edu/~king/


k = length(group);

mat_ind = grp_ind(1) + 1; 
if k > 1
    for ii=2:k 
        mat_ind = mat_ind + grp_ind(ii)*prod(group(1:(ii-1)));
    end
end

end
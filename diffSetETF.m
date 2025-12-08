function [Phi, mat_inds] = diffSetETF(group,inds)
% diffSetETF Create difference set ETF from group elements

% INPUTS:
% group: vector of natural numbers [n_1, ..., n_k] representing a finite 
% % abelian group G = Z_(n_1) x ... x Z_(n_k)
% inds: indices of a difference set in G in an d x k matrix, where d is
% % is the size of the difference set and k is the number of cyclic factors
% % in G. Each ith entry is an integer in the range 0, ..., n_i - 1.

% OUTPUT:
% Phi: synthesis matrix of d x (n_1*n_2*...*n_k) difference set ETF
% mat_inds: the indices in Matlab format of the Fourier matrix used to make Phi

% Author: Emily J King
% https://www.math.colostate.edu/~king/

k = length(group);

% Building lexicographically ordered char table of G
F = fft(eye(group(1)));
if k > 1
    for ii=2:k
        F=kron(F,fft(eye(group(ii))));
    end
end

% Generating Matlab indexing from group indexing
n = prod(group);
d = size(inds,1);
mat_inds = zeros(d,1);
for jj=1:d
     mat_inds(jj) = grp_to_mat(group,inds(jj,:));
end

% Extracting diff set rows from char table to make ETF
Phi = F(mat_inds,:);


end
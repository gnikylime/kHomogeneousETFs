function [PaleyETF,TGal,TCyc,matInds] = genPaleyETF(PaleyDS,p,s)
% genPaleyETF Generate Paley ETF and symmetries

% INPUTS  
% PaleyDS: a (p^s-1)/2 x s int matrix of elements of Z_p x ... x Z_p in the
% % Paley difference set, ordered a^k, k=0, ..., (p^s-3)/2 for a a
% % multiplicative generator of the quadratic residues.
% p: the prime characteristic
% s: the exponent such that the Paley difference set is in F_(p^s)

% OUTPUTS
% PaleyETF: a (p^s-1)/2 x p^s matrix which is the synthesis matrix of a
% % Paley difference set ETF
% TGal: a (p^s-1)/2 x (p^s-1)/2 order s permutation, which maps the rows of
% PaleyETF to their Galois conjugates as quadratic residues in F_(p^s)
% TCyc: a (p^s-1)/2 x (p^s-1)/2 cyclic (order p^s-1)/2) permutation
% matInds: the indices in Matlab format of the Fourier matrix used to make
% % PaleyETF

% Author: Emily J King
% https://www.math.colostate.edu/~king/

if (size(PaleyDS,1) ~= round((p^s-1)/2)) | (size(PaleyDS,2) ~= s) |...
        (mod(p,4) ~= 3)
    error('The inputs do not correspond to a Paley difference set.')
end

%% Generating ETF
d = size(PaleyDS,1);
n = p^s;

[PaleyETF,matInds]=diffSetETF(p*ones(s,1),round(real(PaleyDS)));

%% Generating cyclic translation
TCyc = eye(d);
TCyc = TCyc([end 1:(end-1)], :);

%% Generating Galois permutation

% Row indices of PaleyETF corresponding to quadratic residues in the base
% % field
basefieldQRs=zeros((p-1)/2,1); 

% Sets of s row indices of PaleyETF corresponding to orbits under Galois
% % conjugation
galoisinds=zeros((d-(p-1)/2)/s,s);

% Generate the desired indices
jj = 1; kk = 1;
for ii = 1:d
    tempinds = PaleyDS(ii,:);
    if sum(tempinds(2:end)) == 0
        basefieldQRs(jj) = ii; 
        jj = jj + 1; 
    elseif ~ismember(ii,galoisinds)
        galoisinds(kk,:) = mod((ii-1)*(p.^(0:s-1)),d)+1;
        kk = kk + 1;
    end
end

% Generate Galois permutation 
% Base field quadratic residues are fixed
TGal = zeros(d);
TGal(basefieldQRs,basefieldQRs)=eye((p-1)/2);

% Remaining quadratic residues cycled in groups of s
Ts = eye(s);
Ts = Ts([end 1:end-1], :);
for ii=1:(d-(p-1)/2)/s
    TGal(galoisinds(ii,:),galoisinds(ii,:))=Ts;
end

end



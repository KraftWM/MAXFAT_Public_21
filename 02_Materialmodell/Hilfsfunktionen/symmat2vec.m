function V = symmat2vec(A)
% Speichert relevante Komponenten der symmetrischen "NxN" Matrix A als 
% "N*(N+1)/2 x 1" Vektor
%
%     a11  a12  a13        1  4  6
% A =      a22  a23  ->       2  5 
%               a33              3
%

% ... Größe
N = size(A,1);
% ... Speicher
idx = zeros(N*(N+1)/2,1);
% ... indices
z = 1;
for j = 1 : N
    for i = (j-1):(N-1)
        idx(z) = (i+2-j) + i*N;
        z = z + 1;
    end
end
% ... Vektor
V = A(idx);
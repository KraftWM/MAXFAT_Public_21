function A = vec2symmat(V)
% Speichert relevante Komponenten aus "N*(N+1)/2 x 1" Vektor in
% symmetrische "NxN" Matrix A 
% 
%
%     a11  a12  a13        1  4  6
% A =      a22  a23  ->       2  5 
%               a33              3
%

% ... Größe
N = length(V);
N = (-1+sqrt(1+8*N))/2;
% ... Speicher
A = zeros(N);
idx = zeros(N*(N+1)/2,1);
% ... indices
z = 1;
for j = 1 : N
    for i = (j-1):(N-1)
        idx(z) = (i+2-j) + i*N;
        z = z + 1;
    end
end

% ... Fülle Werte
A(idx) = V;

% ... symmetrischen Teil
for m = 1:N-1
      A(m+1:N, m) = A(m, m+1:N).';
end

end
clear all;

N = 3;

% 1D basis
basis = zeros(N);

for n = 1 : N
     for k = 1 : N
       basis(k, n) = exp(2i * pi * (n-1) * k/N) / sqrt(N);
     end
end


% 2D basis
N1 = 3;
N2 = 3;
basis2D = zeros(N1*N1*N1);

for n = 1 : N1
     for m = 1 : N2
         for l = 1:N
             for k = 1:N
                 basis2D(k, l, m, n) = exp(2i * pi * (k*(n-1)/N1 + l*(m-1)/N2)) / sqrt(N1*N2);
             end
         end
     end
end

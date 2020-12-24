function [index] = stratified(wp)
    N = length(wp);
    index = zeros(1,N);
    Q = cumsum(wp); 
    n = 1;
    m = 1;
    while (n <= N && m <= N)
        u0 = rand/N;
        u = u0 + n/N;
        while (Q(m) < u && m < N)
            m = m+1;
        end
        index(n) = m;
        n = n+1;
    end
end
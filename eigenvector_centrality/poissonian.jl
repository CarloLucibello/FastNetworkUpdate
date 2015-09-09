function first_eig(A; iter = 100)
    v = rand(size(A,1))
    λ=1
    for it=1:iter
        v = A*v
        λ = norm(v)
        v /= λ
    end
    λ, v
end

using LightGraphs
N = 1000
c = 5.
g0 = erdos_renyi(N, c / (N-1), is_directed = false)
A0 = adjacency_matrix(g0)
λ0, v0 = first_eig(A0)
m=zeros(N,N)
for i=1:N
    for j=1:N
        m[i,j] = A0[i,j]
    end
end

Λ0, U0 = eig(m)

g = g0

iter = 1
for it=1:iter
    # EDGE REMOVAL
    i1 = rand(1:N)
    e = out_edges(g, i1)[rand(1:end)]
    rem_edge!(g, e)

    # EDGE ADDITION
    i1, i2 = rand(1:N,2)
    while i1==i2 || has_edge(g, i1, i2)
        i1, i2 = rand(1:N,2)
    end
    add_edge!(g, i1, i2)
end

A = adjacency_matrix(g)
λ, v = first_eig(A)
λ
λ - λ0
norm(v.-v0)
v0*(A-A0)*v0'


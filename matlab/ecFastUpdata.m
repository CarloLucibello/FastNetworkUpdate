c = 5;
nNd = 100;
A = ranSymBinMat(nNd,c/nNd);
pab = [10 6; 6 1]/nNd;
A = blockBinMatrixVec([nNd/2 nNd/2], pab);
g(1:nNd/2) = 1;
g(nNd/2+1:nNd) =2;
A = double(A+A'>0);

spy(A)
[u l] = eig(A);
[l_max pos_max] = max(diag(l));
l0 = l_max;
u0 = u(:,pos_max);
u1 = u0;
eig_est(1) = l_max;
eig_dyna(1) = l_max;

[x y z] = find(A);
L = length(z);
nSteps = 10000;

for s = 1:nSteps
    Aold=A;
    V = zeros(nNd);
    kold = k;
    
    [x y z] = find(A);
    L = length(z);
    k = randi(L);
    rem_i = x(k);
    rem_j = y(k);
    
    if rand() < 1-pab(g(rem_i),g(rem_j))
        removed_link(s,:) = [rem_i rem_j];
        V(x(k),y(k)) = -1;
        V(y(k),x(k)) = -1;
        A(x(k),y(k)) = 0;
        A(y(k),x(k)) = 0;
    end
    
    i = randi(nNd);
    j = randi(nNd);
    count = 1;
    
    while (A(i,j)==1)||(i==j)
        i = randi(nNd);
        j = randi(nNd);
        
        if count > 10^5
            break
        end
        count = count+1;
    end
    
    
    if rand() < pab(g(i),g(j))
        added_link(s,:) = [i j];
        x(k) = i;
        y(k) = j;
        z(k) = 1;
        
        V(i,j) = V(i,j) + 1;
        V(j,i) = V(j,i) + 1;
        A(i,j) = 1;
        A(j,i) = 1;
    end
    
    eig_est(s+1) = eig_est(s) + u1'*V*u1;
    
    [u l] = eig(A);
    [l_max pos_max] = max(diag(l));
    eig_dyna(s+1)= l_max;
    u1 = u(:,pos_max);
    
    distance(s) = min(norm(u0(:)-u1(:)),norm(u0(:)+u1(:)));
    
    density(s) = sum(sum(A))/(nNd*(nNd-1));
    
end
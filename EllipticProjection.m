function [uhI,nodeI,elemI] = EllipticProjection(node,elem,elemVh,uh,h1pGl)
%% nodeI, elemI
nodeI = node(horzcat(elem{:}),:);
elemLen = cellfun('length',elem);
elemI = mat2cell(1:sum(elemLen), 1, elemLen)';

%% uh
chi = cellfun(@(id) uh(id), elemVh, 'UniformOutput', false);    % elementwise d.o.f.s
uhI = cell(size(elem,1),1);

%% uhI
for i = 1:size(elem,1)
    id = elem{i,:};
    X = node(id,:);
    Nv = length(X);
    
    ar_components = X(:,1) .* X([2:Nv,1],2) - X([2:Nv,1],1) .* X(:,2);
    ar = 0.5 * abs(sum(ar_components));                                 % area
    cen = sum((X + X([2:Nv,1],:)) .* repmat(ar_components,1,2))/(6*ar); % centroid
    xK = cen(1,1); yK = cen(1,2);
    hK = 0;                                                             % diameter
    for k = 1:(Nv-1)                                              
        for j = (k+1):Nv                                 
            hK = max(hK, norm(X(k,:) - X(j,:)));    
        end                                                    
    end
    x = node(id,1); y = node(id,2);

    % Basis
    m1 = @(x,y) 1+0*x+0*y;
    m2 = @(x,y) (x-xK)./hK + 0*y; 
    m3 = @(x,y) (y-yK)./hK + 0*x;

    m = {m1, m2, m3};
    Pis = h1pGl{i}; a = Pis*chi{i};
    
    % L2 projection of uh
    uv = zeros(Nv,1);
    for j = 1:3
        uv = uv + a(j)*m{j}(x,y);
    end
    uhI{i} = uv;
end
uhI = cell2mat(uhI);
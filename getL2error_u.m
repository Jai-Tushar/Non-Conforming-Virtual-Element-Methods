function [erruL2] = getL2error_u(node,elem,elemVh,uh,pGl,pde)

chi = cellfun(@(id) uh(id), elemVh, 'UniformOutput', false);    % elementwise d.o.f.s
erruL2 = 0;

for i = 1:size(elem)
    index = elem{i,:};
    X = node(index,:);
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
    
    m1 = @(x) 1+0*x;
    m2 = @(x) (x(1,1)-xK)./hK; 
    m3 = @(x) (x(1,2)-yK)./hK;
    g = {m1, m2, m3};

    
    Pis = pGl{i}; a = Pis*chi{i};
    ag = @(x) 0;  
    for j = 1:3
        ag = @(x) ag(x) + a(j)*g{j}(x);
    end    
    ag   = @(x) ag(x);
    erru = @(x) (pde.exactu(x)-ag(x)).^2;    
    % elementwise error
    nodeT = [node(index,:); cen];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    erruL2 = erruL2 + sum(integralTri(erru,3,nodeT,elemT));

end

erruL2 = sqrt(erruL2);



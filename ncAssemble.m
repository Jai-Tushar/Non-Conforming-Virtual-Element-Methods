%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear non-conforming VEM assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jai Tushar, Postdoctoral fellow, Dept. of Mathematics, IIT Roorkee.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,M,b,h,Freedofs,Dirichlet,PhH1,PhL2] = ncAssemble(order,mesh,pde)

% Dirichlet nodes for a unit square domain 
N = mesh.Pkc4nVh; NT = mesh.Pkn4eVh;
bd = zeros(max(size(N)),min(size(N)));
for i = 1:max(size(N))
    for j = 1:min(size(N))
        if N(i,j) <= 0.1e-8 || N(i,j) >= .99999
            bd(i,j) = i;
        else
            bd(i,j) = 0;
        end
    end
end
Dirichlet = zeros(size(N,1),1);
for i = 1:max(size(N))
    if bd(i,1) == bd(i,2)
        Dirichlet(i,1) = bd(i,1);
    else
        Dirichlet(i,1) = bd(i,1) + bd(i,2);
    end
end
Dirichlet(Dirichlet == 0) = [];

ndof = size(N,1);
nel  = size(NT,1);
fulldofs = 1:ndof;
Freedofs = setdiff(fulldofs,Dirichlet);

A = sparse(ndof,ndof);
M = sparse(ndof,ndof);
b = zeros(ndof,1);
Di = zeros(nel,1);
PhH1 = cell(nel,1);
PhL2 = cell(nel,1);


for el_id = 1:nel
    vert_id = mesh.P1n4e{el_id};
    X = mesh.P1c4n(vert_id,:);
    
    idVh = mesh.Pkn4eVh{el_id};
    
    [ar, ~, hK, G_tilda, G, B, D, H, be] = nclocal(order, X, vert_id, mesh.P1bdEdge, pde);

    h1_p = G\B;
    PhH1{el_id} = h1_p;
   
    h1_st = (eye(length(idVh))- D * h1_p)' * (eye(length(idVh))- D * h1_p);

    Di(el_id,1) = Di(el_id,1) + hK;
    
    % Classical-Stabilization
    AK = (h1_p'*G_tilda*h1_p) + h1_st;

    C = H * (G\B);
    l2_p = H\C;
    PhL2{el_id} = l2_p;
    l2_st = (eye(length(idVh))- D * l2_p)' * (eye(length(idVh))- D * l2_p);
    G_bar = H;
    MK = h1_p'*G_bar*h1_p + ar*l2_st;

    M(idVh,idVh) = M(idVh,idVh) + MK;
    A(idVh,idVh) = A(idVh,idVh) + AK;      

    % Load 
    b(idVh,1) = b(idVh,1) + be; 
    
dbstop if warning
end

h = max(Di);
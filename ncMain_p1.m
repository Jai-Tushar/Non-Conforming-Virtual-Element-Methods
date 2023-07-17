%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear VEM to solve 
%      -\nabla . (K \nabla u) + u = f in \Omega
%                          u  = 0 on \partial \Omega 
% Reference: 1. The nonconforming virtual element method, ESAIM:m2na
%            B. A. de Dios, K. Lipnikov, G. Manzini 
%            2. mVEM package by Y Yu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jai Tushar, Postdoctoral fellow, Dept. of Mathematics, IIT Roorkee.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%% PRE-PROCESSING
% PARAMETERS, DATA and MESH INFO
pde = LEdata;
nr = 5;
kk = 1;
l2err = zeros(nr,1);
h1err = zeros(nr,1);
h = zeros(nr,1);
order = 1;
l2erru = zeros(nr,1);

tic
for k = 1:nr
    % STRUCTURED VORONOI
    if kk == 1
        load(['ncSV_' num2str(k) '.mat']);
        N = mesh.Pkc4nVh; NT = mesh.Pkn4eVh;
    end
    %% PROCESSING
    % Assemble Matrices
    [A,M,b,dia,Freedofs,dirichlet,PhH1,PhL2] = ncAssemble(order, mesh, pde);
    b = b(Freedofs);
    h(k) = dia;
    ndof = size(N,1);
    nel  = size(NT,1);
    
    % Dirichlet boundary conditions
    uD = zeros(ndof,1);
    uD(dirichlet) = pde.exactu(N(dirichlet,:));
    b = b - A(Freedofs,:)*uD;
    % Direct solve
    u = uD;
    Z = (A(Freedofs,Freedofs)+M(Freedofs,Freedofs));
    u(Freedofs) = Z\b;
    %% POST PROCESSING
    [uhI,nodeI,elemI] = EllipticProjection(mesh.P1c4n,mesh.P1n4e,mesh.Pkn4eVh,u,PhH1);
    ueI = pde.exactu(nodeI);

    % Visualise
    figure(2); clf
    set(gcf,'Units','normal');
    set(gcf,'Position',[0.25,0.25,0.6,0.25]);
    subplot(1,2,1)
    plotsol(elemI,nodeI,ueI);
    title('Exact solution')
    xlabel('x'); ylabel('y'); zlabel('u');
    axis tight;
    subplot(1,2,2)
    plotsol(elemI,nodeI,uhI);
    title('Approximate solution')
    xlabel('x'); ylabel('y'); zlabel('u');  

    % Error
    [l2erru(k)] = getL2error_u(mesh.P1c4n,mesh.P1n4e,mesh.Pkn4eVh,u,PhH1,pde);

end
toc

l2order = zeros(nr-1,1);
% h1order = zeros(NR-1,1);
for j=1:nr-1
    l2order(j) = log(l2erru(j)./l2erru(j+1))/log(h(j)/h(j+1));
%     h1order(j) = log(h1err(j)./h1err(j+1))/log(h(j)/h(j+1));
end

fprintf('h = max(diameter)');
h
fprintf('L^2 error and H^1 error');
[l2erru]
fprintf('L^2 EOC and H^1 EOC');
[l2order]

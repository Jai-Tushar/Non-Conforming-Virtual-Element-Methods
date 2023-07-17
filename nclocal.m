% Local Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jai Tushar, Postdoctoral fellow, Dept. of Mathematics, IIT Roorkee.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ar, cen, hK, G_tilda, G, B, D, H, b] = nclocal(order, lc4n, p1n4e, bdEdge, pde)


%-----------Element Geometry-------------
ln4e = {1:length(lc4n)};
Ne = length(lc4n);
ar_components = lc4n(:,1) .* lc4n([2:Ne,1],2) - lc4n([2:Ne,1],1) .* lc4n(:,2);
ar = 0.5 * abs(sum(ar_components));                                       % area
cen = sum((lc4n + lc4n([2:Ne,1],:)) .* repmat(ar_components,1,2))/(6*ar); % centroid
hK = 0;                                                                   % diameter
 for i = 1:(Ne-1)                                              
    for j = (i+1):Ne                                 
        hK = max(hK, norm(lc4n(i,:) - lc4n(j,:)));    
    end                                                    
 end

 clear i j

%-------------Dimensions----------------
nk = 3;
edgeDof = Ne*order;
intDof = order*(order-1)*0.5;
loc_dof = edgeDof + intDof;

%----------------P1-Monomial-Basis-------------------
m1=@(x)(1); m2=@(x)((x(:,1)-cen(1))./hK); m3=@(x)((x(:,2)-cen(2))./hK);
%----------------Gradient-Monomial-Basis-------------
g1=@(x)[0 0]; g2=@(x)[1 0]*(1/hK); g3=@(x)[0 1]*(1/hK); 

%------------Normals------------
temp = circshift(ln4e{:},length(ln4e{:})-1)';
edgeCo = [ln4e{:}' temp];
v1 = 1:Ne; v2 = [2:Ne,1];                                                   % loop index for vertices or edges
[me, wte] = subedgeRoutineGLeg(1,lc4n,edgeCo);   
edgeNo = lc4n(v2,:)-lc4n(v1,:);
he = sqrt(edgeNo(:,1).^2+edgeNo(:,2).^2); 
edgeUnv = [edgeNo(:,2) -edgeNo(:,1)].*(1./he); 

%------------Triangulation----------------------
nodeT = [lc4n(ln4e{:},:); cen];
elemT = [(Ne+1)*ones(Ne,1), (1:Ne)', [2:Ne,1]'];

%------------ G_tilda Matrix -------------------
intm = @(x) [dot(g1(x),g1(x)) dot(g1(x),g2(x)) dot(g1(x),g3(x));
             dot(g2(x),g1(x)) dot(g2(x),g2(x)) dot(g2(x),g3(x));
             dot(g3(x),g1(x)) dot(g3(x),g2(x)) dot(g3(x),g3(x))];

G_tilda = integralTri(intm,1,nodeT,elemT);      

clear intm

%------------ D Matrix -------------------
D = zeros(loc_dof,nk);

mod_wrap = @(x,a)mod(x-1,a) + 1;
m = @(x) [m1(x) m2(x) m3(x)];

% Gauss Legendre Quadrature
for i = 1:edgeDof
    ma = m(me(mod_wrap(i,Ne),:));  
    D(i,:) = D(i,:) + 0.5*wte*(ma);
end

clear i m ma

%------------------B Matrix----------------
B_tilda = zeros(nk,loc_dof);

gm = @(x) [g1(x) g2(x) g3(x)];


for i = 1:edgeDof
    gma = reshape(gm(me(mod_wrap(i,Ne),:)),2,[])';   
    ne = edgeUnv(mod_wrap(i,Ne),:);
    gmadotne = sum(gma.*ne,2);
    B_tilda(:,i) = B_tilda(:,i) + he(i)*gmadotne;
    
    clear gmadotne j
end

B = B_tilda;
B(1,:) = he;

clear i j gm

%-----------------G Matrix----------------
G = B*D;

%-----------------H Matrix----------------
m = @(x) [m1(x)*m1(x) m1(x)*m2(x) m1(x)*m3(x);
          m2(x)*m1(x) m2(x)*m2(x) m2(x)*m3(x);
          m3(x)*m1(x) m3(x)*m2(x) m3(x)*m3(x)];

H = integralTri(m,2,nodeT,elemT);    


%------------ Load --------------
f = pde.f;

intf = integralTri(f,3,nodeT,elemT);
b = (1/loc_dof)*intf*ones(loc_dof,1);
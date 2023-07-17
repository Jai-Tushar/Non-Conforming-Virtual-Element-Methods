%% Basic Principles / Equiv Projector
% function b_local = loadterm(X,pde)
%     n_sides = length(X);
%     [Area,centroid] = localGBDH(X);
%     b_local = (1/n_sides)*Area*(pde.f(centroid));
% end

% % % function b_local = loadterm(X,f)
% % %     n_sides = length(X);
% % %     [Area,centroid] = localGBDH(X);
% % %     b_local = (1/n_sides)*Area*(f(centroid));
% % % end

%% Hitchikers Guide
function b = loadterm(X,f)

dof = length(X);
[~,cen,hd,G,B,~,H]=localGBDH(X);
h1_proj = (G\eye(size(G)))*B;
C = H * (G\eye(size(G))) * B;
l2_proj = (H\eye(size(H))) * C;
wt = TriGaussPoints(3);

int_1 = 0; int_2 = 0; int_3 = 0;
for i = 1:dof
    Xx = [X(i,1) X(mod(i,dof)+1,1) cen(1)];
    Yy = [X(i,2) X(mod(i,dof)+1,2) cen(2)];
    
    z1 = 0; z2 = 0; z3 = 0;
    for j = 1:length(wt(:,1))
        xn = Xx(1)*(1-wt(j,1)-wt(j,2)) + Xx(2)*wt(j,1) + Xx(3)*wt(j,2);
        yn = Yy(1)*(1-wt(j,1)-wt(j,2)) + Yy(2)*wt(j,1) + Yy(3)*wt(j,2);
        
        v = [xn,yn];
        z1 = z1 + f(v)*wt(j,3);
        z2 = z2 + f(v)*((xn - cen(1))/hd)*wt(j,3);
        z3 = z3 + f(v)*((yn - cen(2))/hd)*wt(j,3);
    end
    int_1 = int_1 + polyarea(Xx,Yy)*z1;
    int_2 = int_2 + polyarea(Xx,Yy)*z2;
    int_3 = int_3 + polyarea(Xx,Yy)*z3;
end
% % b = l2_proj(1,:)*int_1 + l2_proj(2,:)*int_2 + l2_proj(3,:)*int_3;
% % b = h1_proj(1,:)*int_1 + h1_proj(2,:)*int_2 + h1_proj(3,:)*int_3;
b = [int_1 int_2 int_3]*l2_proj;
end
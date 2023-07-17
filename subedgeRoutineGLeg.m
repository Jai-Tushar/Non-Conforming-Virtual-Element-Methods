% Edge Routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jai Tushar, PhD, Dept. of Mathematics, BITS-Goa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edgeDofMat, Wt] = subedgeRoutineGLeg(order,lc4n,edgeCo)

Ne = length(lc4n);
edgeDof = Ne*(order);
edgeDofMat = zeros(edgeDof,2);

[xGL, Wt] = GaussLegendre(order);
xGL = flip(xGL);


for j = 1:Ne
    xi = lc4n(edgeCo(j,1),1); yi = lc4n(edgeCo(j,1),2);
    xj = lc4n(edgeCo(j,2),1); yj = lc4n(edgeCo(j,2),2);
    pt = zeros(order,2);
    for i = 1:length(xGL)
        t = 0.5*xGL(i) + 0.5;
        pt(i,1) = xi + t*(xj-xi);
        pt(i,2) = yi + t*(yj-yi);
    end
    for k = 1:order
        edgeDofMat((k-1)*Ne+j,:) = pt(k,:);
    end
end
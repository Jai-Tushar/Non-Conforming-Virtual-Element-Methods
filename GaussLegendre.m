% Gauss Legendre qudrature: Weights and Nodes

function [x, W] = GaussLegendre(n)

% Legendre polynomial
p=polegende(n);

% Polynomial roots
x=roots(p(n+1,:));

% Polynomial derivative
pn=polyder(p(n+1,:));

% Weights
for i=1:n
   W(i)=2./((1-x(i).^2).*((polyval(pn,x(i))).^2));
end
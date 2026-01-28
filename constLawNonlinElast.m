function [ W, dWdG, dWdu, d2Wd2G, d2Wd2u, d2WdGdu, d3Wd3G, d3Wd2Gdu, d3WdGd2u ] = constLawNonlinElast( G, u, elK, elG, trf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dudx = G(1:2);
dudy = G(3:4);

grad_u = [ dudx.'  dudy.' ];
grad_u_e = zeros(3);
grad_u_e(1:2,1:2) = grad_u;

iFp = diag( [ 1/trf 1/trf 1 ] );

I = eye(3);
F = I + grad_u_e;
Fe = F * iFp;
Be = Fe * Fe.';
Ce = Fe.' * Fe;
Ee = ( Ce - I )/2;
Se = (elK-2*elG/3)*I*tensTrace3(Ee) + 2*elG*Ee;

iFpT4I = tens2mlt4I(iFp.');
Fe4I = tens2mlt4I(Fe);
FeT4I = tens2mlt4I(Fe.');
Be4I = tens2mlt4I(Be);
Se4I = tens2mlt4I(Se);
vFe = tens2vec(Fe);

C3 = zeros(9);
C3(1:3,1:3) = eye(3);
C3(4:6,7:9) = eye(3);
C3(7:9,4:6) = eye(3);

W = (1/2)*tensTrace3(Se*Ee);
dWdFe = Fe*Se;
d2Wd2Fe = C3*Se4I + (elK-2*elG/3)*(vFe*vFe.') + elG*Be4I*C3 + elG*(Fe4I*FeT4I*C3);

dWdF = dWdFe*(iFp.');
vdWdF = tens2vec(dWdF);
dWdG = vdWdF([1 7 4 2]).';
d2Wd2F = iFpT4I.' * C3 * d2Wd2Fe * C3 * iFpT4I;
d2Wd2G_m = d2Wd2F([1 7 4 2],[1 7 4 2]);
d2Wd2G = d2Wd2G_m(:).';

d3Wd3G = zeros(1,64);
d3Wd3G(1) = F(1,1)*(3*elK+4*elG)*iFp(1,1)^4;
d3Wd3G(22) = F(2,1)*(3*elK+4*elG)*iFp(1,1)^4;
d3Wd3G(43) = F(1,2)*(3*elK+4*elG)*iFp(2,2)^4;
d3Wd3G(64) = F(2,2)*(3*elK+4*elG)*iFp(2,2)^4;
d3Wd3G([21 18  6]) = F(1,1)*(elK+4*elG/3)*iFp(1,1)^4;
d3Wd3G([41 35 11]) = F(1,1)*(elK+4*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([17  5  2]) = F(2,1)*(elK+4*elG/3)*iFp(1,1)^4;
d3Wd3G([62 56 32]) = F(2,1)*(elK+4*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([33  9  3]) = F(1,2)*(elK+4*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([63 60 48]) = F(1,2)*(elK+4*elG/3)*iFp(2,2)^4;
d3Wd3G([24 30 54]) = F(2,2)*(elK+4*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([44 47 59]) = F(2,2)*(elK+4*elG/3)*iFp(2,2)^4;
d3Wd3G([58 46 55 31 28 40]) = F(1,1)*elG*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([57 45 51 15 12 36]) = F(2,1)*elG*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([53 29 50 14  8 20]) = F(1,2)*elG*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([10 34  7 19 25 37]) = F(2,2)*elG*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([61 52 16]) = F(1,1)*(elK-2*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([42 39 27]) = F(2,1)*(elK-2*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([38 23 26]) = F(1,2)*(elK-2*elG/3)*iFp(1,1)^2*iFp(2,2)^2;
d3Wd3G([49  4 13]) = F(2,2)*(elK-2*elG/3)*iFp(1,1)^2*iFp(2,2)^2;

dWdu = zeros(1,2);
d2Wd2u = zeros(1,4);
d2WdGdu = zeros(1,8);
d3Wd2Gdu = zeros(1,32);
d3WdGd2u = zeros(1,16);

end


function [ W, dWdG, dWdu, d2Wd2G, d2Wd2u, d2WdGdu, d3Wd3G, d3Wd2Gdu, d3WdGd2u ] = constLawMagDem( G, u, mgA, mgK, ani, cf1, cf2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dmdx = G(1:3);
dmdy = G(6:8);

dldx = G(4);
dldy = G(9);

dndx = G(5);
dndy = G(10);

m = u(1:3);
lam = u(4);

%% magnetism

p = [ cos(ani) sin(ani) 0 ];

W_mag = (1/2)*mgA*( dmdx*dmdx.' + dmdy*dmdy.' ) - (1/2)*mgK*(p*m.')^2 + lam*( m*m.' - 1 );

%% field

grad_n_e = [ dndx dndy 0 ];

W_field = (-1/2)*cf1*( grad_n_e * grad_n_e.' ) + cf2*( m * grad_n_e.' );

regPar = 0;
W_reg = (1/2)*regPar*( dldx^2 + dldy^2 );

W = W_mag + W_field + W_reg;

%% deriv. first

dWdG = zeros(1,10);
dWdG(1,1:3) = mgA*dmdx;
dWdG(1,6:8) = mgA*dmdy;
dWdG(1,5) = -cf1*dndx + cf2*m(1);
dWdG(1,10) = -cf1*dndy + cf2*m(2);
dWdG(1,4) = regPar*dldx;
dWdG(1,9) = regPar*dldy;

dWdu = zeros(1,5);
dWdu(1,1:3) = -mgK*(p*m.')*p + 2*lam*m + cf2*grad_n_e;
dWdu(1,4) = m*m.' - 1;

%% deriv. second

d2Wd2G_r = zeros(10,10);
d2Wd2G_r(1:3,1:3) = mgA*eye(3);
d2Wd2G_r(6:8,6:8) = mgA*eye(3);
d2Wd2G_r(5,5) = -cf1;
d2Wd2G_r(10,10) = -cf1;
d2Wd2G_r(4,4) = regPar;
d2Wd2G_r(9,9) = regPar;
d2Wd2G = d2Wd2G_r(:).';

d2Wd2u_r = zeros(5,5);
d2Wd2u_r(1:3,1:3) = -mgK*(p.')*p + 2*lam*eye(3);
d2Wd2u_r(4,1:3) = 2*m;
d2Wd2u_r(1:3,4) = 2*m.';
d2Wd2u = d2Wd2u_r(:).';

d2WdGdu = zeros(1,50);
d2WdGdu(1,5) = cf2;
d2WdGdu(1,20) = cf2;

%% deriv. third

d3Wd3G = zeros(1,1000);
d3Wd2Gdu = zeros(1,500);
d3WdGd2u = zeros(1,250);

end


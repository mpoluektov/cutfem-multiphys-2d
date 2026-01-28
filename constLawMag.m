function [ W, dWdG, dWdu, d2Wd2G, d2Wd2u, d2WdGdu, d3Wd3G, d3Wd2Gdu, d3WdGd2u ] = constLawMag( G, u, mgA, mgK, ani )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dmdx = G(1:3);
dmdy = G(5:7);
m = u(1:3);
lam = u(4);

p = [ cos(ani) sin(ani) 0 ];

W = (1/2)*mgA*( dmdx*dmdx.' + dmdy*dmdy.' ) - (1/2)*mgK*(p*m.')^2 + lam*( m*m.' - 1 );

%% deriv. first

dWdG = zeros(1,8);
dWdG(1,1:3) = mgA*dmdx;
dWdG(1,5:7) = mgA*dmdy;

dWdu = zeros(1,4);
dWdu(1,1:3) = -mgK*(p*m.')*p + 2*lam*m;
dWdu(1,4) = m*m.' - 1;

%% deriv. second

d2Wd2G_r = zeros(8,8);
d2Wd2G_r(1:3,1:3) = mgA*eye(3);
d2Wd2G_r(5:7,5:7) = mgA*eye(3);
d2Wd2G = d2Wd2G_r(:).';

d2Wd2u_r = zeros(4,4);
d2Wd2u_r(1:3,1:3) = -mgK*(p.')*p + 2*lam*eye(3);
d2Wd2u_r(4,1:3) = 2*m;
d2Wd2u_r(1:3,4) = 2*m.';
d2Wd2u = d2Wd2u_r(:).';

d2WdGdu = zeros(1,32);

%% deriv. third

d3Wd3G = zeros(1,512);
d3Wd2Gdu = zeros(1,256);
d3WdGd2u = zeros(1,128);

end


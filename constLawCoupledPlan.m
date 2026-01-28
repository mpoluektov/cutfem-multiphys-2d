function [ W, dWdG, dWdu, d2Wd2G, d2Wd2u, d2WdGdu, d3Wd3G, d3Wd2Gdu, d3WdGd2u ] = constLawCoupledPlan( G, u, elK, elG, trf, mgA, mgK, ani, cf1, cf2, cf3, SDLT )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dudx = G(1:2);
dudy = G(5:6);

dtdx = G(3);
dtdy = G(7);

dndx = G(4);
dndy = G(8);

t = u(3);
m = [ cos(t) sin(t) 0 ];
dmdt = [ -sin(t) cos(t) 0 ];
d2md2t = [ -cos(t) -sin(t) 0 ];

%% mechanics

grad_u = [ dudx.' dudy.' ];
grad_u_e = zeros(3);
grad_u_e(1:2,1:2) = grad_u;

iFp = diag( [ 1/trf trf 1 ] );

I = eye(3);
F = I + grad_u_e;
J = det(F);
Fe = F * iFp;
Be = Fe * Fe.';
Ce = Fe.' * Fe;
Ee = ( Ce - I )/2;
Se = (elK-2*elG/3)*I*tensTrace3(Ee) + 2*elG*Ee;

W_mech = (1/2)*tensTrace3(Se*Ee);

%% magnetism

p = [ cos(ani) sin(ani) 0 ];

iF = tensInv3(F);
W_mag = (1/2)*mgA*( dtdx^2 + dtdy^2 ) - (1/2)*mgK*(m*iF.'*p.')^2;

%% field

grad_n_e = [ dndx dndy 0 ].';

iC = iF * iF.';
W_field = (-1/2) * cf1 * J * grad_n_e.' * iC * grad_n_e + cf2 * m * iF.' * grad_n_e + (-1/2) * cf3 / J;

W = W_mech + W_mag + W_field;

%% deriv. first

FT4I = tens2mlt4I(F.');
iFpT4I = tens2mlt4I(iFp.');
iF4I = tens2mlt4I(iF);
iFT4I = tens2mlt4I(iF.');
iC4I = tens2mlt4I(iC);
Fe4I = tens2mlt4I(Fe);
FeT4I = tens2mlt4I(Fe.');
Be4I = tens2mlt4I(Be);
Se4I = tens2mlt4I(Se);
viFT = tens2vec(iF.');
vFe = tens2vec(Fe);

mp = ( iF.' * p.' ) * ( iF * m.' ).';
mp_der = ( iF.' * p.' ) * ( iF * dmdt.' ).';
vmp = tens2vec(mp);
mpmp = vmp * vmp.';
mp4I = tens2mlt4I(mp);
mpT4I = tens2mlt4I(mp.');

C1 = eye(9);
C3 = zeros(9);
C3(1:3,1:3) = eye(3);
C3(4:6,7:9) = eye(3);
C3(7:9,4:6) = eye(3);

der1 = cf1 * J * iC * ( grad_n_e * grad_n_e.' ) * iF - (1/2) * cf1 * J * iF * ( grad_n_e.' * iC * grad_n_e ) - cf2 * iF * ( m.' * grad_n_e.' ) * iF + (1/2) * cf3 * iF / J;
der2 = -cf1 * J * iC * grad_n_e + cf2 * iF * m.';

dWdF = Fe*Se*(iFp.') + mgK*(m*iF.'*p.')*mp + der1.';
vdWdF = tens2vec(dWdF);

dWdG = zeros(1,8);
dWdG(1,1:2) = vdWdF([1 7]).';
dWdG(1,5:6) = vdWdF([4 2]).';
dWdG(1,3) = mgA*dtdx;
dWdG(1,7) = mgA*dtdy;
dWdG(1,4) = der2(1);
dWdG(1,8) = der2(2);

dWdu = zeros(1,4);
dWdu(1,3) = -mgK*(m*iF.'*p.')*(dmdt*iF.'*p.') + cf2 * dmdt * iF.' * grad_n_e;

%% deriv. second

d2Wd2Fe = C3*Se4I + (elK-2*elG/3)*(vFe*vFe.') + elG*Be4I*C3 + elG*(Fe4I*FeT4I*C3);

interimTens1 = iC * ( grad_n_e * grad_n_e.' ) * iF;
interimTens2 = ( grad_n_e * grad_n_e.' ) * iF;
der3 = -cf1 * J * C3 * tens2mlt4I(interimTens1) * iFT4I + cf1 * J * tens2vec(interimTens1.') * viFT.' - ...
    cf1 * J * tens2mlt4I(interimTens2.') * C3 * iC4I * iC4I * (C1+C3) * FT4I * C3 + ...
    cf1 * J * viFT * tens2vec(interimTens1.').' + (-1/2) * cf1 * J * ( grad_n_e.' * iC * grad_n_e ) * ( viFT * viFT.' ) + ...
    (1/2) * cf1 * J * ( grad_n_e.' * iC * grad_n_e ) * C3 * iF4I * iFT4I + ...
    cf2 * C3 * tens2mlt4I( iF * ( m.' * grad_n_e.' ) ) * C3 * iF4I * iFT4I + cf2 * tens2mlt4I( iF.' * ( m.' * grad_n_e.' ).' ) * iF4I * iFT4I + ...
    (-1/2) * cf3 * ( viFT * viFT.' ) / J + (-1/2) * cf3 * C3 * iF4I * iFT4I / J;

d2Wd2F = iFpT4I.' * C3 * d2Wd2Fe * C3 * iFpT4I - mgK*mpmp - mgK*(m*iF.'*p.')*(mp4I*iF4I*C3) - mgK*(m*iF.'*p.')*(iFT4I*mpT4I*C3) + der3;

vtrm1 = -cf2 * iF * m.' + cf1 * J * iC * grad_n_e;
vtrm2 = -cf1 * J * iC * grad_n_e;
vtrm3 = cf1 * J * iF.' * grad_n_e;

d2Wd2G_r = zeros(8,8);
d2Wd2G_r([1 2 5 6],[1 2 5 6]) = d2Wd2F([1 7 4 2],[1 7 4 2]);
d2Wd2G_r(3,3) = mgA;
d2Wd2G_r(7,7) = mgA;
d2Wd2G_r(4,4) = -cf1 * J * iC(1,1);
d2Wd2G_r(4,8) = -cf1 * J * iC(1,2);
d2Wd2G_r(8,4) = -cf1 * J * iC(2,1);
d2Wd2G_r(8,8) = -cf1 * J * iC(2,2);
d2Wd2G_r(4,1) = vtrm1(1) * iF(1,1) + iF(1,1) * vtrm2(1) + iC(1,1) * vtrm3(1);
d2Wd2G_r(4,2) = vtrm1(1) * iF(1,2) + iF(1,2) * vtrm2(1) + iC(1,1) * vtrm3(2);
d2Wd2G_r(1,4) = vtrm1(1) * iF(1,1) + iF(1,1) * vtrm2(1) + iC(1,1) * vtrm3(1);
d2Wd2G_r(2,4) = vtrm1(1) * iF(1,2) + iF(1,2) * vtrm2(1) + iC(1,1) * vtrm3(2);
d2Wd2G_r(8,1) = vtrm1(1) * iF(2,1) + iF(1,1) * vtrm2(2) + iC(1,2) * vtrm3(1);
d2Wd2G_r(8,2) = vtrm1(1) * iF(2,2) + iF(1,2) * vtrm2(2) + iC(1,2) * vtrm3(2);
d2Wd2G_r(1,8) = vtrm1(1) * iF(2,1) + iF(1,1) * vtrm2(2) + iC(1,2) * vtrm3(1);
d2Wd2G_r(2,8) = vtrm1(1) * iF(2,2) + iF(1,2) * vtrm2(2) + iC(1,2) * vtrm3(2);
d2Wd2G_r(4,5) = vtrm1(2) * iF(1,1) + iF(2,1) * vtrm2(1) + iC(1,2) * vtrm3(1);
d2Wd2G_r(4,6) = vtrm1(2) * iF(1,2) + iF(2,2) * vtrm2(1) + iC(1,2) * vtrm3(2);
d2Wd2G_r(5,4) = vtrm1(2) * iF(1,1) + iF(2,1) * vtrm2(1) + iC(1,2) * vtrm3(1);
d2Wd2G_r(6,4) = vtrm1(2) * iF(1,2) + iF(2,2) * vtrm2(1) + iC(1,2) * vtrm3(2);
d2Wd2G_r(8,5) = vtrm1(2) * iF(2,1) + iF(2,1) * vtrm2(2) + iC(2,2) * vtrm3(1);
d2Wd2G_r(8,6) = vtrm1(2) * iF(2,2) + iF(2,2) * vtrm2(2) + iC(2,2) * vtrm3(2);
d2Wd2G_r(5,8) = vtrm1(2) * iF(2,1) + iF(2,1) * vtrm2(2) + iC(2,2) * vtrm3(1);
d2Wd2G_r(6,8) = vtrm1(2) * iF(2,2) + iF(2,2) * vtrm2(2) + iC(2,2) * vtrm3(2);
d2Wd2G = d2Wd2G_r(:).';

d2Wd2u_r = zeros(4,4);
d2Wd2u_r(3,3) = -mgK*(m*iF.'*p.')*(d2md2t*iF.'*p.') - mgK*(dmdt*iF.'*p.')*(dmdt*iF.'*p.') + cf2 * d2md2t * iF.' * grad_n_e;
d2Wd2u = d2Wd2u_r(:).';

der4 = -cf2 * iF * ( dmdt.' * grad_n_e.' ) * iF;
der5 = cf2 * iF * dmdt.';

dWdF_der = mgK*(m*iF.'*p.')*mp_der + mgK*(dmdt*iF.'*p.')*mp + der4.';
vdWdF_der = tens2vec(dWdF_der);

d2WdGdu_r = zeros(8,4);
d2WdGdu_r([1 2],3) = vdWdF_der([1 7]).';
d2WdGdu_r([5 6],3) = vdWdF_der([4 2]).';
d2WdGdu_r(4,3) = der5(1);
d2WdGdu_r(8,3) = der5(2);
d2WdGdu = d2WdGdu_r(:).';

%% deriv. third

d3Wd3G = zeros(1,8*8*8);
d3Wd2Gdu = zeros(1,8*8*4);
d3WdGd2u = zeros(1,8*4*4);

if ( nargout == 9 )
    d3Wd3G_num = zeros(8*8,8);
    d3Wd2Gdu_num = zeros(8*8,4);
    d3WdGd2u_num = zeros(8*4,4);
    for ii=1:8
        G_c = G;
        G_c(ii) = G_c(ii) + SDLT;
        [ dum1, dum2, dum3, d2Wd2G_c, dum4, dum5 ] = constLawCoupledPlan( G_c, u, elK, elG, trf, mgA, mgK, ani, cf1, cf2, cf3 );
        d3Wd3G_num(:,ii) = ( d2Wd2G_c - d2Wd2G ).' / SDLT;
    end
    d3Wd3G = d3Wd3G_num(:).';
    for ii=1:4
        u_c = u;
        u_c(ii) = u_c(ii) + SDLT;
        [ dum1, dum2, dum3, d2Wd2G_c, dum4, d2WdGdu_c ] = constLawCoupledPlan( G, u_c, elK, elG, trf, mgA, mgK, ani, cf1, cf2, cf3 );
        d3Wd2Gdu_num(:,ii) = ( d2Wd2G_c - d2Wd2G ).' / SDLT;
        d3WdGd2u_num(:,ii) = ( d2WdGdu_c - d2WdGdu ).' / SDLT;
    end
    d3Wd2Gdu = d3Wd2Gdu_num(:).';
    d3WdGd2u = d3WdGd2u_num(:).';
end

end


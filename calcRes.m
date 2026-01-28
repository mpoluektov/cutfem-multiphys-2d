function [ E, F, J, resElem ] = calcRes( u, matPar, geomPar, numPar, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

elKU = matPar.elKU;
elKT = matPar.elKT;
elGU = matPar.elGU;
elGT = matPar.elGT;
trfU = matPar.trfU;
trfT = matPar.trfT;
mgAU = matPar.mgAU;
mgAT = matPar.mgAT;
mgKU = matPar.mgKU;
mgKT = matPar.mgKT;
aniU = matPar.aniU;
aniT = matPar.aniT;
cf1U = matPar.cf1U;
cf1T = matPar.cf1T;
cf2U = matPar.cf2U;
cf2T = matPar.cf2T;
cf3U = matPar.cf3U;
cf3T = matPar.cf3T;

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;
dofs = geomPar.dofs;

SDLT = numPar.SDLT;
gam = numPar.NIT_PAR;
kap = numPar.STAB_PAR;

u_U = u( 1:(dofs*Nx*Ny), : );
u_T = u( (dofs*Nx*Ny+1):(2*dofs*Nx*Ny), : );

E = 0;

F_U = zeros( Nx*Ny, dofs );
F_T = zeros( Nx*Ny, dofs );

J_UU = zeros( dofs*Nx*Ny, dofs*Nx*Ny );
J_TT = zeros( dofs*Nx*Ny, dofs*Nx*Ny );
J_UT = zeros( dofs*Nx*Ny, dofs*Nx*Ny );
J_TU = zeros( dofs*Nx*Ny, dofs*Nx*Ny );

%% get interface elements and nodes

intElem = intParam.intElem;
intConn = intParam.intConn;
cutLen = intParam.cutLen;
elemNorms = intParam.elemNorms;
fracElem = intParam.fracElem;
fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;
elTypes = intParam.elTypes;
exclNodes_U = intParam.exclNodes_U;
exclNodes_T = intParam.exclNodes_T;
intNodesUnk_C = intParam.intNodesUnk_C;

%% calculate energies and their derivatives in elements

W_U_elem_od = zeros( (Nx-1)*(Ny-1), 4 );
W_U_elem_ev = zeros( (Nx-1)*(Ny-1), 4 );
W_T_elem_od = zeros( (Nx-1)*(Ny-1), 4 );
W_T_elem_ev = zeros( (Nx-1)*(Ny-1), 4 );

G_U_elem_od = zeros( (Nx-1)*(Ny-1), 2*dofs );
G_U_elem_ev = zeros( (Nx-1)*(Ny-1), 2*dofs );
G_T_elem_od = zeros( (Nx-1)*(Ny-1), 2*dofs );
G_T_elem_ev = zeros( (Nx-1)*(Ny-1), 2*dofs );

u_U_elem_od = zeros( (Nx-1)*(Ny-1), dofs );
u_U_elem_ev = zeros( (Nx-1)*(Ny-1), dofs );
u_T_elem_od = zeros( (Nx-1)*(Ny-1), dofs );
u_T_elem_ev = zeros( (Nx-1)*(Ny-1), dofs );

dWdG_U_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs );
dWdG_U_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs );
dWdG_T_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs );
dWdG_T_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs );

dWdu_U_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs );
dWdu_U_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs );
dWdu_T_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs );
dWdu_T_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs );

d2Wd2G_U_elem_od = zeros( (Nx-1)*(Ny-1), 16*dofs*dofs );
d2Wd2G_U_elem_ev = zeros( (Nx-1)*(Ny-1), 16*dofs*dofs );
d2Wd2G_T_elem_od = zeros( (Nx-1)*(Ny-1), 16*dofs*dofs );
d2Wd2G_T_elem_ev = zeros( (Nx-1)*(Ny-1), 16*dofs*dofs );

d2Wd2u_U_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs );
d2Wd2u_U_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs );
d2Wd2u_T_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs );
d2Wd2u_T_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs );

d2WdGdu_U_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs );
d2WdGdu_U_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs );
d2WdGdu_T_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs );
d2WdGdu_T_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs );

d3Wd3G_U_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs*dofs );
d3Wd3G_U_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs*dofs );
d3Wd3G_T_elem_od = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs*dofs );
d3Wd3G_T_elem_ev = zeros( (Nx-1)*(Ny-1), 8*dofs*dofs*dofs );

d3Wd2Gdu_U_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs*dofs );
d3Wd2Gdu_U_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs*dofs );
d3Wd2Gdu_T_elem_od = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs*dofs );
d3Wd2Gdu_T_elem_ev = zeros( (Nx-1)*(Ny-1), 4*dofs*dofs*dofs );

d3WdGd2u_U_elem_od = zeros( (Nx-1)*(Ny-1), 2*dofs*dofs*dofs );
d3WdGd2u_U_elem_ev = zeros( (Nx-1)*(Ny-1), 2*dofs*dofs*dofs );
d3WdGd2u_T_elem_od = zeros( (Nx-1)*(Ny-1), 2*dofs*dofs*dofs );
d3WdGd2u_T_elem_ev = zeros( (Nx-1)*(Ny-1), 2*dofs*dofs*dofs );

elemNodeInds_elem_od = zeros( (Nx-1)*(Ny-1), 4 );
elemNodeInds_elem_ev = zeros( (Nx-1)*(Ny-1), 4 );

sfts = (dofs-1):(-1):0;

for ii=1:((Nx-1)*(Ny-1))

    ind_bl = ii + floor((ii-0.5)/(Nx-1));
    ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
    ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
    ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));
    
    elem = [ ind_tl  ind_bl  ind_br   1 ;
             ind_br  ind_tr  ind_tl  -1 ];
    
    for jj=1:2
        
        elInd = 2*ii+jj-2;
        
        if ( elTypes(elInd) == 2 )||( elTypes(elInd) == 3 )
            %. untransformed or intersected

            Gx = elem(jj,4) * ( u_U(dofs*elem(jj,3)-sfts) - u_U(dofs*elem(jj,2)-sfts) ).' * (1/dh);
            Gy = elem(jj,4) * ( u_U(dofs*elem(jj,1)-sfts) - u_U(dofs*elem(jj,2)-sfts) ).' * (1/dh);
            G_loc = [ Gx Gy ];
            u_locA = ( u_U(dofs*elem(jj,1)-sfts) + u_U(dofs*elem(jj,2)-sfts) + u_U(dofs*elem(jj,3)-sfts) ).' * (1/3);
            u_locB = ( 3*u_U(dofs*elem(jj,1)-sfts) + u_U(dofs*elem(jj,2)-sfts) + u_U(dofs*elem(jj,3)-sfts) ).' * (1/5);
            u_locC = ( u_U(dofs*elem(jj,1)-sfts) + 3*u_U(dofs*elem(jj,2)-sfts) + u_U(dofs*elem(jj,3)-sfts) ).' * (1/5);
            u_locD = ( u_U(dofs*elem(jj,1)-sfts) + u_U(dofs*elem(jj,2)-sfts) + 3*u_U(dofs*elem(jj,3)-sfts) ).' * (1/5);

            constLaw = 1;
            if ( constLaw == 1 )
                %. nonlinear elastic constitutive law
                [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawNonlinElast( G_loc, u_locA, elKU, elGU, trfU );
                W_locB = W_locA;
                W_locC = W_locA;
                W_locD = W_locA;
                dWdG_locB = dWdG_locA;
                dWdG_locC = dWdG_locA;
                dWdG_locD = dWdG_locA;
                dWdu_locB = dWdu_locA;
                dWdu_locC = dWdu_locA;
                dWdu_locD = dWdu_locA;
                d2Wd2G_locB = d2Wd2G_locA;
                d2Wd2G_locC = d2Wd2G_locA;
                d2Wd2G_locD = d2Wd2G_locA;
                d2Wd2u_locB = d2Wd2u_locA;
                d2Wd2u_locC = d2Wd2u_locA;
                d2Wd2u_locD = d2Wd2u_locA;
                d2WdGdu_locB = d2WdGdu_locA;
                d2WdGdu_locC = d2WdGdu_locA;
                d2WdGdu_locD = d2WdGdu_locA;
            elseif ( constLaw == 2 )
                %. coupled magneto-mechanics constitutive law
                if ( elTypes(elInd) == 3 )
                    [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawCoupledPlan( G_loc, u_locA, elKU, elGU, trfU, mgAU, mgKU, aniU, cf1U, cf2U, cf3U, SDLT );
                else
                    [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA ] = constLawCoupledPlan( G_loc, u_locA, elKU, elGU, trfU, mgAU, mgKU, aniU, cf1U, cf2U, cf3U, SDLT );
                end
                [ W_locB, dWdG_locB, dWdu_locB, d2Wd2G_locB, d2Wd2u_locB, d2WdGdu_locB ] = constLawCoupledPlan( G_loc, u_locB, elKU, elGU, trfU, mgAU, mgKU, aniU, cf1U, cf2U, cf3U, SDLT );
                [ W_locC, dWdG_locC, dWdu_locC, d2Wd2G_locC, d2Wd2u_locC, d2WdGdu_locC ] = constLawCoupledPlan( G_loc, u_locC, elKU, elGU, trfU, mgAU, mgKU, aniU, cf1U, cf2U, cf3U, SDLT );
                [ W_locD, dWdG_locD, dWdu_locD, d2Wd2G_locD, d2Wd2u_locD, d2WdGdu_locD ] = constLawCoupledPlan( G_loc, u_locD, elKU, elGU, trfU, mgAU, mgKU, aniU, cf1U, cf2U, cf3U, SDLT );
            elseif ( constLaw == 3 )
                %. micromagnetism constitutive law
                [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawMag( G_loc, u_locA, mgAU, mgKU, aniU );
                [ W_locB, dWdG_locB, dWdu_locB, d2Wd2G_locB, d2Wd2u_locB, d2WdGdu_locB ] = constLawMag( G_loc, u_locB, mgAU, mgKU, aniU );
                [ W_locC, dWdG_locC, dWdu_locC, d2Wd2G_locC, d2Wd2u_locC, d2WdGdu_locC ] = constLawMag( G_loc, u_locC, mgAU, mgKU, aniU );
                [ W_locD, dWdG_locD, dWdu_locD, d2Wd2G_locD, d2Wd2u_locD, d2WdGdu_locD ] = constLawMag( G_loc, u_locD, mgAU, mgKU, aniU );
            end

            if ( jj == 1 )
                W_U_elem_od(ii,:) = [ W_locA W_locB W_locC W_locD ];
                G_U_elem_od(ii,:) = G_loc;
                u_U_elem_od(ii,:) = u_locA;
                dWdG_U_elem_od(ii,:) = [ dWdG_locA dWdG_locB dWdG_locC dWdG_locD ];
                dWdu_U_elem_od(ii,:) = [ dWdu_locA dWdu_locB dWdu_locC dWdu_locD ];
                d2Wd2G_U_elem_od(ii,:) = [ d2Wd2G_locA d2Wd2G_locB d2Wd2G_locC d2Wd2G_locD ];
                d2Wd2u_U_elem_od(ii,:) = [ d2Wd2u_locA d2Wd2u_locB d2Wd2u_locC d2Wd2u_locD ];
                d2WdGdu_U_elem_od(ii,:) = [ d2WdGdu_locA d2WdGdu_locB d2WdGdu_locC d2WdGdu_locD ];
                if ( elTypes(elInd) == 3 )
                    d3Wd3G_U_elem_od(ii,:) = d3Wd3G_locA;
                    d3Wd2Gdu_U_elem_od(ii,:) = d3Wd2Gdu_locA;
                    d3WdGd2u_U_elem_od(ii,:) = d3WdGd2u_locA;
                end
            else
                W_U_elem_ev(ii,:) = [ W_locA W_locB W_locC W_locD ];
                G_U_elem_ev(ii,:) = G_loc;
                u_U_elem_ev(ii,:) = u_locA;
                dWdG_U_elem_ev(ii,:) = [ dWdG_locA dWdG_locB dWdG_locC dWdG_locD ];
                dWdu_U_elem_ev(ii,:) = [ dWdu_locA dWdu_locB dWdu_locC dWdu_locD ];
                d2Wd2G_U_elem_ev(ii,:) = [ d2Wd2G_locA d2Wd2G_locB d2Wd2G_locC d2Wd2G_locD ];
                d2Wd2u_U_elem_ev(ii,:) = [ d2Wd2u_locA d2Wd2u_locB d2Wd2u_locC d2Wd2u_locD ];
                d2WdGdu_U_elem_ev(ii,:) = [ d2WdGdu_locA d2WdGdu_locB d2WdGdu_locC d2WdGdu_locD ];
                if ( elTypes(elInd) == 3 )
                    d3Wd3G_U_elem_ev(ii,:) = d3Wd3G_locA;
                    d3Wd2Gdu_U_elem_ev(ii,:) = d3Wd2Gdu_locA;
                    d3WdGd2u_U_elem_ev(ii,:) = d3WdGd2u_locA;
                end
            end
        end
        if ( elTypes(elInd) == 1 )||( elTypes(elInd) == 3 )
            %. transformed or intersected

            Gx = elem(jj,4) * ( u_T(dofs*elem(jj,3)-sfts) - u_T(dofs*elem(jj,2)-sfts) ).' * (1/dh);
            Gy = elem(jj,4) * ( u_T(dofs*elem(jj,1)-sfts) - u_T(dofs*elem(jj,2)-sfts) ).' * (1/dh);
            G_loc = [ Gx Gy ];
            u_locA = ( u_T(dofs*elem(jj,1)-sfts) + u_T(dofs*elem(jj,2)-sfts) + u_T(dofs*elem(jj,3)-sfts) ).' * (1/3);
            u_locB = ( 3*u_T(dofs*elem(jj,1)-sfts) + u_T(dofs*elem(jj,2)-sfts) + u_T(dofs*elem(jj,3)-sfts) ).' * (1/5);
            u_locC = ( u_T(dofs*elem(jj,1)-sfts) + 3*u_T(dofs*elem(jj,2)-sfts) + u_T(dofs*elem(jj,3)-sfts) ).' * (1/5);
            u_locD = ( u_T(dofs*elem(jj,1)-sfts) + u_T(dofs*elem(jj,2)-sfts) + 3*u_T(dofs*elem(jj,3)-sfts) ).' * (1/5);

            constLaw = 1;
            if ( constLaw == 1 )
                %. nonlinear elastic constitutive law
                [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawNonlinElast( G_loc, u_locA, elKT, elGT, trfT );
                W_locB = W_locA;
                W_locC = W_locA;
                W_locD = W_locA;
                dWdG_locB = dWdG_locA;
                dWdG_locC = dWdG_locA;
                dWdG_locD = dWdG_locA;
                dWdu_locB = dWdu_locA;
                dWdu_locC = dWdu_locA;
                dWdu_locD = dWdu_locA;
                d2Wd2G_locB = d2Wd2G_locA;
                d2Wd2G_locC = d2Wd2G_locA;
                d2Wd2G_locD = d2Wd2G_locA;
                d2Wd2u_locB = d2Wd2u_locA;
                d2Wd2u_locC = d2Wd2u_locA;
                d2Wd2u_locD = d2Wd2u_locA;
                d2WdGdu_locB = d2WdGdu_locA;
                d2WdGdu_locC = d2WdGdu_locA;
                d2WdGdu_locD = d2WdGdu_locA;
            elseif ( constLaw == 2 )
                %. coupled magneto-mechanics constitutive law
                if ( elTypes(elInd) == 3 )
                    [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawCoupledPlan( G_loc, u_locA, elKT, elGT, trfT, mgAT, mgKT, aniT, cf1T, cf2T, cf3T, SDLT );
                else
                    [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA ] = constLawCoupledPlan( G_loc, u_locA, elKT, elGT, trfT, mgAT, mgKT, aniT, cf1T, cf2T, cf3T, SDLT );
                end
                [ W_locB, dWdG_locB, dWdu_locB, d2Wd2G_locB, d2Wd2u_locB, d2WdGdu_locB ] = constLawCoupledPlan( G_loc, u_locB, elKT, elGT, trfT, mgAT, mgKT, aniT, cf1T, cf2T, cf3T, SDLT );
                [ W_locC, dWdG_locC, dWdu_locC, d2Wd2G_locC, d2Wd2u_locC, d2WdGdu_locC ] = constLawCoupledPlan( G_loc, u_locC, elKT, elGT, trfT, mgAT, mgKT, aniT, cf1T, cf2T, cf3T, SDLT );
                [ W_locD, dWdG_locD, dWdu_locD, d2Wd2G_locD, d2Wd2u_locD, d2WdGdu_locD ] = constLawCoupledPlan( G_loc, u_locD, elKT, elGT, trfT, mgAT, mgKT, aniT, cf1T, cf2T, cf3T, SDLT );
            elseif ( constLaw == 3 )
                %. micromagnetism constitutive law
                [ W_locA, dWdG_locA, dWdu_locA, d2Wd2G_locA, d2Wd2u_locA, d2WdGdu_locA, d3Wd3G_locA, d3Wd2Gdu_locA, d3WdGd2u_locA ] = constLawMag( G_loc, u_locA, mgAT, mgKT, aniT );
                [ W_locB, dWdG_locB, dWdu_locB, d2Wd2G_locB, d2Wd2u_locB, d2WdGdu_locB ] = constLawMag( G_loc, u_locB, mgAT, mgKT, aniT );
                [ W_locC, dWdG_locC, dWdu_locC, d2Wd2G_locC, d2Wd2u_locC, d2WdGdu_locC ] = constLawMag( G_loc, u_locC, mgAT, mgKT, aniT );
                [ W_locD, dWdG_locD, dWdu_locD, d2Wd2G_locD, d2Wd2u_locD, d2WdGdu_locD ] = constLawMag( G_loc, u_locD, mgAT, mgKT, aniT );
            end

            if ( jj == 1 )
                W_T_elem_od(ii,:) = [ W_locA W_locB W_locC W_locD ];
                G_T_elem_od(ii,:) = G_loc;
                u_T_elem_od(ii,:) = u_locA;
                dWdG_T_elem_od(ii,:) = [ dWdG_locA dWdG_locB dWdG_locC dWdG_locD ];
                dWdu_T_elem_od(ii,:) = [ dWdu_locA dWdu_locB dWdu_locC dWdu_locD ];
                d2Wd2G_T_elem_od(ii,:) = [ d2Wd2G_locA d2Wd2G_locB d2Wd2G_locC d2Wd2G_locD ];
                d2Wd2u_T_elem_od(ii,:) = [ d2Wd2u_locA d2Wd2u_locB d2Wd2u_locC d2Wd2u_locD ];
                d2WdGdu_T_elem_od(ii,:) = [ d2WdGdu_locA d2WdGdu_locB d2WdGdu_locC d2WdGdu_locD ];
                if ( elTypes(elInd) == 3 )
                    d3Wd3G_T_elem_od(ii,:) = d3Wd3G_locA;
                    d3Wd2Gdu_T_elem_od(ii,:) = d3Wd2Gdu_locA;
                    d3WdGd2u_T_elem_od(ii,:) = d3WdGd2u_locA;
                end
            else
                W_T_elem_ev(ii,:) = [ W_locA W_locB W_locC W_locD ];
                G_T_elem_ev(ii,:) = G_loc;
                u_T_elem_ev(ii,:) = u_locA;
                dWdG_T_elem_ev(ii,:) = [ dWdG_locA dWdG_locB dWdG_locC dWdG_locD ];
                dWdu_T_elem_ev(ii,:) = [ dWdu_locA dWdu_locB dWdu_locC dWdu_locD ];
                d2Wd2G_T_elem_ev(ii,:) = [ d2Wd2G_locA d2Wd2G_locB d2Wd2G_locC d2Wd2G_locD ];
                d2Wd2u_T_elem_ev(ii,:) = [ d2Wd2u_locA d2Wd2u_locB d2Wd2u_locC d2Wd2u_locD ];
                d2WdGdu_T_elem_ev(ii,:) = [ d2WdGdu_locA d2WdGdu_locB d2WdGdu_locC d2WdGdu_locD ];
                if ( elTypes(elInd) == 3 )
                    d3Wd3G_T_elem_ev(ii,:) = d3Wd3G_locA;
                    d3Wd2Gdu_T_elem_ev(ii,:) = d3Wd2Gdu_locA;
                    d3WdGd2u_T_elem_ev(ii,:) = d3WdGd2u_locA;
                end
            end
        end
        if ( jj == 1 )
            elemNodeInds_elem_od(ii,:) = elem(jj,:);
        else
            elemNodeInds_elem_ev(ii,:) = elem(jj,:);
        end
    end
end

W_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = W_U_elem_ev;
W_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = W_U_elem_od;
W_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = W_T_elem_ev;
W_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = W_T_elem_od;

G_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = G_U_elem_ev;
G_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = G_U_elem_od;
G_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = G_T_elem_ev;
G_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = G_T_elem_od;

u_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = u_U_elem_ev;
u_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = u_U_elem_od;
u_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = u_T_elem_ev;
u_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = u_T_elem_od;

dWdG_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = dWdG_U_elem_ev;
dWdG_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = dWdG_U_elem_od;
dWdG_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = dWdG_T_elem_ev;
dWdG_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = dWdG_T_elem_od;

dWdu_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = dWdu_U_elem_ev;
dWdu_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = dWdu_U_elem_od;
dWdu_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = dWdu_T_elem_ev;
dWdu_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = dWdu_T_elem_od;

d2Wd2G_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2Wd2G_U_elem_ev;
d2Wd2G_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2Wd2G_U_elem_od;
d2Wd2G_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2Wd2G_T_elem_ev;
d2Wd2G_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2Wd2G_T_elem_od;

d2Wd2u_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2Wd2u_U_elem_ev;
d2Wd2u_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2Wd2u_U_elem_od;
d2Wd2u_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2Wd2u_T_elem_ev;
d2Wd2u_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2Wd2u_T_elem_od;

d2WdGdu_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2WdGdu_U_elem_ev;
d2WdGdu_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2WdGdu_U_elem_od;
d2WdGdu_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d2WdGdu_T_elem_ev;
d2WdGdu_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d2WdGdu_T_elem_od;

d3Wd3G_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3Wd3G_U_elem_ev;
d3Wd3G_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3Wd3G_U_elem_od;
d3Wd3G_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3Wd3G_T_elem_ev;
d3Wd3G_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3Wd3G_T_elem_od;

d3Wd2Gdu_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3Wd2Gdu_U_elem_ev;
d3Wd2Gdu_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3Wd2Gdu_U_elem_od;
d3Wd2Gdu_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3Wd2Gdu_T_elem_ev;
d3Wd2Gdu_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3Wd2Gdu_T_elem_od;

d3WdGd2u_U_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3WdGd2u_U_elem_ev;
d3WdGd2u_U_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3WdGd2u_U_elem_od;
d3WdGd2u_T_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d3WdGd2u_T_elem_ev;
d3WdGd2u_T_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d3WdGd2u_T_elem_od;

elemNodeInds_elem( 2:2:(2*(Nx-1)*(Ny-1)), : ) = elemNodeInds_elem_ev;
elemNodeInds_elem( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = elemNodeInds_elem_od;

%% export

resElem.W_U_elem = W_U_elem;
resElem.W_T_elem = W_T_elem;
resElem.G_U_elem = G_U_elem;
resElem.G_T_elem = G_T_elem;
resElem.u_U_elem = u_U_elem;
resElem.u_T_elem = u_T_elem;
resElem.dWdG_U_elem = dWdG_U_elem;
resElem.dWdG_T_elem = dWdG_T_elem;
resElem.dWdu_U_elem = dWdu_U_elem;
resElem.dWdu_T_elem = dWdu_T_elem;
resElem.d2Wd2G_U_elem = d2Wd2G_U_elem;
resElem.d2Wd2G_T_elem = d2Wd2G_T_elem;
resElem.d2Wd2u_U_elem = d2Wd2u_U_elem;
resElem.d2Wd2u_T_elem = d2Wd2u_T_elem;
resElem.d2WdGdu_U_elem = d2WdGdu_U_elem;
resElem.d2WdGdu_T_elem = d2WdGdu_T_elem;

%% calculate integrals

gradv = [ 1  0 ;
         -1 -1 ;
          0  1 ;
          0 -1 ;
         -1  0 ;
          1  1 ];
dGdu_ref = [ zeros(dofs)  -eye(dofs)  eye(dofs)   ;
             eye(dofs)    -eye(dofs)  zeros(dofs) ];
elem_edeges = [ 1  7  1 ;
                1  6  2 ;
                4  1  3 ;
                8  4  2 ;
                2  4  1 ;
                9  2  3 ;
                2  5  2 ;
               10  5  1 ;
                5  3  3 ;
                3 11  2 ;
                3  6  1 ;
                6 12  3 ];
dudu_ref = [ (1/3)*eye(dofs)  (1/3)*eye(dofs)  (1/3)*eye(dofs) ;
             (3/5)*eye(dofs)  (1/5)*eye(dofs)  (1/5)*eye(dofs) ;
             (1/5)*eye(dofs)  (3/5)*eye(dofs)  (1/5)*eye(dofs) ;
             (1/5)*eye(dofs)  (1/5)*eye(dofs)  (3/5)*eye(dofs) ];
normals = [ 1          0         ;
            0          1         ;
            1/sqrt(2)  1/sqrt(2) ];
lengths = [ 1  1  sqrt(2) ];
sz = [ dofs*Nx*Ny  dofs*Nx*Ny ];
valv = [ 1/3 1/5 1/5 3/5 ;
         1/3 1/5 3/5 1/5 ;
         1/3 3/5 1/5 1/5 ;
         1/3 3/5 1/5 1/5 ;
         1/3 1/5 1/5 3/5 ;
         1/3 1/5 3/5 1/5 ];
pweight = [ -27 25 25 25 ]/48;

for ii=1:(Nx*Ny)
               
    %. squares
    sq_ar = ii - floor( (ii-0.5)/Nx );
    sq_al = ii - floor( (ii-0.5)/Nx ) - 1;
    sq_br = ii - floor( (ii-0.5)/Nx ) - (Nx-1);
    sq_bl = ii - floor( (ii-0.5)/Nx ) - (Nx-1) - 1;

    sq_all = ii - floor( (ii-0.5)/Nx ) - 2;
    sq_aal = ii - floor( (ii-0.5)/Nx ) + (Nx-1) - 1;
    sq_brr = ii - floor( (ii-0.5)/Nx ) - (Nx-1) + 1;
    sq_bbr = ii - floor( (ii-0.5)/Nx ) - 2*(Nx-1);
    
    elem_inds = [ 2*sq_al-1 ;
                  2*sq_ar-1 ;
                  2*sq_br-1 ;
                    2*sq_al ;
                    2*sq_br ;
                    2*sq_bl ;
                   2*sq_all ;
                 2*sq_aal-1 ;
                    2*sq_ar ;
                 2*sq_brr-1 ;
                   2*sq_bbr ;
                  2*sq_bl-1 ];

    if ( ii == 1 )
        %. LB corner
        inda = 2;
    elseif ( ii == Nx*Ny-Nx+1 )
        %. LT corner
        inda = [ 3 5 ];
    elseif ( ii == Nx )
        %. RB corner
        inda = [ 1 4 ];
    elseif ( ii == Nx*Ny )
        %. RT corner
        inda = 6;
    elseif ( rem(ii,Nx) == 1 )
        %. L side
        inda = [ 2 3 5 ];
    elseif ( ii > Nx*Ny-Nx+1 )
        %. T side
        inda = [ 3 5 6 ];
    elseif ( rem(ii,Nx) == 0 )
        %. R side
        inda = [ 1 4 6 ];
    elseif ( ii < Nx )
        %. B side
        inda = [ 1 2 4 ];
    else
        %. inside
        inda = [ 1 2 3 4 5 6 ];
    end
    
    %. find edges
    indmsk = ones(1,12);
    if ( rem(ii,Nx) == 1 )
        %. L side
        indmsk([1 2 3 4 5 11 12]) = 0;
    end
    if ( rem(ii,Nx) == 2 )
        %. L side dist. by h
        indmsk(1) = 0;
    end
    if ( rem(ii,Nx) == 0 )
        %. R side
        indmsk([5 6 7 8 9 10 11]) = 0;
    end
    if ( rem(ii,Nx) == Nx-1 )
        %. R side dist. by h
        indmsk(8) = 0;
    end
    if ( ii <= Nx )
        %. B side
        indmsk([2 7 8 9 10 11 12]) = 0;
    end
    if ( ii <= 2*Nx )&&( ii >= Nx+1 )
        %. B side dist. by h
        indmsk(10) = 0;
    end
    if ( ii >= Nx*Ny-Nx+1 )
        %. T side
        indmsk([1 2 3 4 5 6 7]) = 0;
    end
    if ( ii >= Nx*Ny-2*Nx+1 )&&( ii <= Nx*Ny-Nx )
        %. T side dist. by h
        indmsk(4) = 0;
    end
    inde = find( indmsk == 1 );
    
    for jj=1:size(inda,2)
        
        dWdu_U = reshape( dWdu_U_elem( elem_inds(inda(jj)), : ), dofs, 4 );
        dWdu_T = reshape( dWdu_T_elem( elem_inds(inda(jj)), : ), dofs, 4 );
        dWdG_U = reshape( dWdG_U_elem( elem_inds(inda(jj)), : ), dofs, 8 );
        dWdG_T = reshape( dWdG_T_elem( elem_inds(inda(jj)), : ), dofs, 8 );
        d2Wd2G_U = reshape( d2Wd2G_U_elem( elem_inds(inda(jj)), : ), 2*dofs, 8*dofs );
        d2Wd2G_T = reshape( d2Wd2G_T_elem( elem_inds(inda(jj)), : ), 2*dofs, 8*dofs );
        d2Wd2u_U = reshape( d2Wd2u_U_elem( elem_inds(inda(jj)), : ), dofs, 4*dofs );
        d2Wd2u_T = reshape( d2Wd2u_T_elem( elem_inds(inda(jj)), : ), dofs, 4*dofs );
        d2WdGdu_U = reshape( d2WdGdu_U_elem( elem_inds(inda(jj)), : ), 2*dofs, 4*dofs );
        d2WdGdu_T = reshape( d2WdGdu_T_elem( elem_inds(inda(jj)), : ), 2*dofs, 4*dofs );
        d3Wd3G_U = reshape( d3Wd3G_U_elem( elem_inds(inda(jj)), : ), 4*dofs*dofs, 2*dofs );
        d3Wd3G_T = reshape( d3Wd3G_T_elem( elem_inds(inda(jj)), : ), 4*dofs*dofs, 2*dofs );
        d3Wd2Gdu_U = reshape( d3Wd2Gdu_U_elem( elem_inds(inda(jj)), : ), 4*dofs*dofs, dofs );
        d3Wd2Gdu_T = reshape( d3Wd2Gdu_T_elem( elem_inds(inda(jj)), : ), 4*dofs*dofs, dofs );
        d3WdGd2u_U = reshape( d3WdGd2u_U_elem( elem_inds(inda(jj)), : ), 2*dofs*dofs, dofs );
        d3WdGd2u_T = reshape( d3WdGd2u_T_elem( elem_inds(inda(jj)), : ), 2*dofs*dofs, dofs );

        elem = elemNodeInds_elem( elem_inds(inda(jj)), : );
        elType_loc = elTypes( elem_inds(inda(jj)) );
        if ( elType_loc == 2 )
            %. untransformed
            area_U = 1;
            area_T = 0;
        elseif ( elType_loc == 1 )
            %. transformed
            area_U = 0;
            area_T = 1;
        else
            %. intersected
            indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
            area_U = fracElem(indEl_loc);
            area_T = 1 - fracElem(indEl_loc);
        end
        
        gphi = gradv(inda(jj),:).' * (1/dh);
        cphi = valv(inda(jj),:).' .* pweight.';
        int_U = area_U * ( dWdG_U * kron(pweight,eye(2)).' * gphi + dWdu_U * cphi ) * (1/2) * (dh^2);
        int_T = area_T * ( dWdG_T * kron(pweight,eye(2)).' * gphi + dWdu_T * cphi ) * (1/2) * (dh^2);
        F_U(ii,:) = F_U(ii,:) + int_U.';
        F_T(ii,:) = F_T(ii,:) + int_T.';
        
        %. N.R. derivatives
        rw = dofs*ii-sfts;
        rows = repmat(rw.',1,3*dofs);
        idx = [ dofs*elem(1)-sfts  dofs*elem(2)-sfts  dofs*elem(3)-sfts ];
        cols = repmat(idx,dofs,1);
        ind = sub2ind(sz,rows,cols);

        dGdu = elem(4) * dGdu_ref * (1/dh);
        dudu_cut = dudu_ref(1:dofs,:);

        pweight_e = reshape(repmat(pweight,dofs,1),4*dofs,1);
        der_int_U = area_U * kron(gphi.',eye(dofs)) * ( d2Wd2G_U * kron(pweight,eye(2*dofs)).' * dGdu + d2WdGdu_U * diag(pweight_e) * dudu_ref ) * (1/2) * (dh^2);
        der_int_T = area_T * kron(gphi.',eye(dofs)) * ( d2Wd2G_T * kron(pweight,eye(2*dofs)).' * dGdu + d2WdGdu_T * diag(pweight_e) * dudu_ref ) * (1/2) * (dh^2);
        J_UU(ind(:)) = J_UU(ind(:)) + der_int_U(:);
        J_TT(ind(:)) = J_TT(ind(:)) + der_int_T(:);

        cphi_e = reshape(repmat(cphi.',dofs,1),4*dofs,1);
        der_int_U = area_U * ( kron(cphi.',eye(dofs)) * d2WdGdu_U.' * dGdu + d2Wd2u_U * diag(cphi_e) * dudu_ref ) * (1/2) * (dh^2);
        der_int_T = area_T * ( kron(cphi.',eye(dofs)) * d2WdGdu_T.' * dGdu + d2Wd2u_T * diag(cphi_e) * dudu_ref ) * (1/2) * (dh^2);
        J_UU(ind(:)) = J_UU(ind(:)) + der_int_U(:);
        J_TT(ind(:)) = J_TT(ind(:)) + der_int_T(:);
        
        %. line integrals
        if ( elType_loc == 3 )

            %. intersected
            indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
            elIntPts = intConn(indEl_loc,:);
            elInt_p1n1ind = intNodes( elIntPts(1), 1 );
            elInt_p1n1w = fracFaces( elIntPts(1) );
            elInt_p1n2ind = intNodes( elIntPts(1), 2 );
            elInt_p1n2w = 1 - fracFaces( elIntPts(1) );
            elInt_p2n1ind = intNodes( elIntPts(2), 1 );
            elInt_p2n1w = fracFaces( elIntPts(2) );
            elInt_p2n2ind = intNodes( elIntPts(2), 2 );
            elInt_p2n2w = 1 - fracFaces( elIntPts(2) );

            phi_p1n1 = 0;
            if ( elInt_p1n1ind == ii )
                phi_p1n1 = 1;
            end
            phi_p1n2 = 0;
            if ( elInt_p1n2ind == ii )
                phi_p1n2 = 1;
            end
            phi_p2n1 = 0;
            if ( elInt_p2n1ind == ii )
                phi_p2n1 = 1;
            end
            phi_p2n2 = 0;
            if ( elInt_p2n2ind == ii )
                phi_p2n2 = 1;
            end
            phi_p1 = phi_p1n1 * elInt_p1n1w + phi_p1n2 * elInt_p1n2w;
            phi_p2 = phi_p2n1 * elInt_p2n1w + phi_p2n2 * elInt_p2n2w;

            elInt_p1n1idx = dofs*elInt_p1n1ind-sfts;
            elInt_p1n2idx = dofs*elInt_p1n2ind-sfts;
            elInt_p2n1idx = dofs*elInt_p2n1ind-sfts;
            elInt_p2n2idx = dofs*elInt_p2n2ind-sfts;

            uU_p1n1 = u_U( elInt_p1n1idx );
            uU_p1n2 = u_U( elInt_p1n2idx );
            uU_p2n1 = u_U( elInt_p2n1idx );
            uU_p2n2 = u_U( elInt_p2n2idx );
            uU_p1 = uU_p1n1 * elInt_p1n1w + uU_p1n2 * elInt_p1n2w;
            uU_p2 = uU_p2n1 * elInt_p2n1w + uU_p2n2 * elInt_p2n2w;

            uT_p1n1 = u_T( elInt_p1n1idx );
            uT_p1n2 = u_T( elInt_p1n2idx );
            uT_p2n1 = u_T( elInt_p2n1idx );
            uT_p2n2 = u_T( elInt_p2n2idx );
            uT_p1 = uT_p1n1 * elInt_p1n1w + uT_p1n2 * elInt_p1n2w;
            uT_p2 = uT_p2n1 * elInt_p2n1w + uT_p2n2 * elInt_p2n2w;

            uJ_p1 = ( uT_p1 - uU_p1 );
            uJ_p2 = ( uT_p2 - uU_p2 );

            secLen = cutLen( indEl_loc );
            secNorm = elemNorms( indEl_loc, : ).';

            %. first integral: [phi].dW.n
            int_A = secLen * ( phi_p1 + phi_p2 ) * ( dWdG_U(:,1:2) + dWdG_T(:,1:2) ) * secNorm * (1/4);
            F_U(ii,:) = F_U(ii,:) + int_A.';
            F_T(ii,:) = F_T(ii,:) - int_A.';

            %. first integral, deriv.
            der_int_AU = secLen * ( phi_p1 + phi_p2 ) * kron(secNorm.',eye(dofs)) * ( d2Wd2G_U(:,1:(2*dofs)) * dGdu + d2WdGdu_U(:,1:dofs) * dudu_cut ) * (1/4);
            der_int_AT = secLen * ( phi_p1 + phi_p2 ) * kron(secNorm.',eye(dofs)) * ( d2Wd2G_T(:,1:(2*dofs)) * dGdu + d2WdGdu_T(:,1:dofs) * dudu_cut ) * (1/4);
            J_UU(ind(:)) = J_UU(ind(:)) + der_int_AU(:);
            J_TT(ind(:)) = J_TT(ind(:)) - der_int_AT(:);
            J_UT(ind(:)) = J_UT(ind(:)) + der_int_AT(:);
            J_TU(ind(:)) = J_TU(ind(:)) - der_int_AU(:);
            
            %. second integral part one: n[u]:d2W:dphi
            uJ_norm = ( uJ_p1 + uJ_p2 ) * secNorm.';
            uJ_norm_v = reshape( uJ_norm, 2*dofs, 1 );
            der_U_v = d2Wd2G_U(:,1:(2*dofs)) * uJ_norm_v;
            der_T_v = d2Wd2G_T(:,1:(2*dofs)) * uJ_norm_v;
            der_U = reshape( der_U_v, dofs, 2 );
            der_T = reshape( der_T_v, dofs, 2 );
            int_U = secLen * der_U * gphi * (1/4);
            int_T = secLen * der_T * gphi * (1/4);
            F_U(ii,:) = F_U(ii,:) - int_U.';
            F_T(ii,:) = F_T(ii,:) - int_T.';

            %. second integral part two: n[u]:d2W:dphi
            vphi = valv(inda(jj),1);
            der_U = d2WdGdu_U(:,1:dofs).' * uJ_norm_v;
            der_T = d2WdGdu_T(:,1:dofs).' * uJ_norm_v;
            int_U = secLen * der_U * vphi * (1/4);
            int_T = secLen * der_T * vphi * (1/4);
            F_U(ii,:) = F_U(ii,:) - int_U.';
            F_T(ii,:) = F_T(ii,:) - int_T.';

            %. second integral part one, deriv.
            der_int_U = secLen * kron(gphi.',eye(dofs)) * kron(eye(2*dofs),uJ_norm_v.') * ( d3Wd3G_U * dGdu + d3Wd2Gdu_U * dudu_cut ) * (1/4);
            der_int_T = secLen * kron(gphi.',eye(dofs)) * kron(eye(2*dofs),uJ_norm_v.') * ( d3Wd3G_T * dGdu + d3Wd2Gdu_T * dudu_cut ) * (1/4);
            der_int_U_p1n1w = kron(gphi.',eye(dofs)) * d2Wd2G_U(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p1n1w * secLen * (1/4);
            der_int_U_p1n2w = kron(gphi.',eye(dofs)) * d2Wd2G_U(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p1n2w * secLen * (1/4);
            der_int_U_p2n1w = kron(gphi.',eye(dofs)) * d2Wd2G_U(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p2n1w * secLen * (1/4);
            der_int_U_p2n2w = kron(gphi.',eye(dofs)) * d2Wd2G_U(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p2n2w * secLen * (1/4);
            der_int_T_p1n1w = kron(gphi.',eye(dofs)) * d2Wd2G_T(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p1n1w * secLen * (1/4);
            der_int_T_p1n2w = kron(gphi.',eye(dofs)) * d2Wd2G_T(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p1n2w * secLen * (1/4);
            der_int_T_p2n1w = kron(gphi.',eye(dofs)) * d2Wd2G_T(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p2n1w * secLen * (1/4);
            der_int_T_p2n2w = kron(gphi.',eye(dofs)) * d2Wd2G_T(:,1:(2*dofs)) * kron(secNorm,eye(dofs)) * elInt_p2n2w * secLen * (1/4);
            J_UU(ind(:)) = J_UU(ind(:)) - der_int_U(:);
            J_TT(ind(:)) = J_TT(ind(:)) - der_int_T(:);
            J_UU(rw,elInt_p1n1idx) = J_UU(rw,elInt_p1n1idx) + der_int_U_p1n1w;
            J_UU(rw,elInt_p1n2idx) = J_UU(rw,elInt_p1n2idx) + der_int_U_p1n2w;
            J_UU(rw,elInt_p2n1idx) = J_UU(rw,elInt_p2n1idx) + der_int_U_p2n1w;
            J_UU(rw,elInt_p2n2idx) = J_UU(rw,elInt_p2n2idx) + der_int_U_p2n2w;
            J_UT(rw,elInt_p1n1idx) = J_UT(rw,elInt_p1n1idx) - der_int_U_p1n1w;
            J_UT(rw,elInt_p1n2idx) = J_UT(rw,elInt_p1n2idx) - der_int_U_p1n2w;
            J_UT(rw,elInt_p2n1idx) = J_UT(rw,elInt_p2n1idx) - der_int_U_p2n1w;
            J_UT(rw,elInt_p2n2idx) = J_UT(rw,elInt_p2n2idx) - der_int_U_p2n2w;
            J_TT(rw,elInt_p1n1idx) = J_TT(rw,elInt_p1n1idx) - der_int_T_p1n1w;
            J_TT(rw,elInt_p1n2idx) = J_TT(rw,elInt_p1n2idx) - der_int_T_p1n2w;
            J_TT(rw,elInt_p2n1idx) = J_TT(rw,elInt_p2n1idx) - der_int_T_p2n1w;
            J_TT(rw,elInt_p2n2idx) = J_TT(rw,elInt_p2n2idx) - der_int_T_p2n2w;
            J_TU(rw,elInt_p1n1idx) = J_TU(rw,elInt_p1n1idx) + der_int_T_p1n1w;
            J_TU(rw,elInt_p1n2idx) = J_TU(rw,elInt_p1n2idx) + der_int_T_p1n2w;
            J_TU(rw,elInt_p2n1idx) = J_TU(rw,elInt_p2n1idx) + der_int_T_p2n1w;
            J_TU(rw,elInt_p2n2idx) = J_TU(rw,elInt_p2n2idx) + der_int_T_p2n2w;

            %. second integral part two, deriv.
            der_int_U_a = secLen * vphi * kron(eye(dofs),uJ_norm_v.') * ( d3WdGd2u_U * dudu_cut ) * (1/4);
            der_int_T_a = secLen * vphi * kron(eye(dofs),uJ_norm_v.') * ( d3WdGd2u_T * dudu_cut ) * (1/4);
            der_int_U_b = secLen * vphi * d3Wd2Gdu_U.' * kron(eye(2*dofs),uJ_norm_v) * dGdu * (1/4);
            der_int_T_b = secLen * vphi * d3Wd2Gdu_T.' * kron(eye(2*dofs),uJ_norm_v) * dGdu * (1/4);
            der_int_U_p1n1w = vphi * d2WdGdu_U(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p1n1w * secLen * (1/4);
            der_int_U_p1n2w = vphi * d2WdGdu_U(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p1n2w * secLen * (1/4);
            der_int_U_p2n1w = vphi * d2WdGdu_U(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p2n1w * secLen * (1/4);
            der_int_U_p2n2w = vphi * d2WdGdu_U(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p2n2w * secLen * (1/4);
            der_int_T_p1n1w = vphi * d2WdGdu_T(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p1n1w * secLen * (1/4);
            der_int_T_p1n2w = vphi * d2WdGdu_T(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p1n2w * secLen * (1/4);
            der_int_T_p2n1w = vphi * d2WdGdu_T(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p2n1w * secLen * (1/4);
            der_int_T_p2n2w = vphi * d2WdGdu_T(:,1:dofs).' * kron(secNorm,eye(dofs)) * elInt_p2n2w * secLen * (1/4);
            J_UU(ind(:)) = J_UU(ind(:)) - der_int_U_a(:) - der_int_U_b(:);
            J_TT(ind(:)) = J_TT(ind(:)) - der_int_T_a(:) - der_int_T_b(:);
            J_UU(rw,elInt_p1n1idx) = J_UU(rw,elInt_p1n1idx) + der_int_U_p1n1w;
            J_UU(rw,elInt_p1n2idx) = J_UU(rw,elInt_p1n2idx) + der_int_U_p1n2w;
            J_UU(rw,elInt_p2n1idx) = J_UU(rw,elInt_p2n1idx) + der_int_U_p2n1w;
            J_UU(rw,elInt_p2n2idx) = J_UU(rw,elInt_p2n2idx) + der_int_U_p2n2w;
            J_UT(rw,elInt_p1n1idx) = J_UT(rw,elInt_p1n1idx) - der_int_U_p1n1w;
            J_UT(rw,elInt_p1n2idx) = J_UT(rw,elInt_p1n2idx) - der_int_U_p1n2w;
            J_UT(rw,elInt_p2n1idx) = J_UT(rw,elInt_p2n1idx) - der_int_U_p2n1w;
            J_UT(rw,elInt_p2n2idx) = J_UT(rw,elInt_p2n2idx) - der_int_U_p2n2w;
            J_TT(rw,elInt_p1n1idx) = J_TT(rw,elInt_p1n1idx) - der_int_T_p1n1w;
            J_TT(rw,elInt_p1n2idx) = J_TT(rw,elInt_p1n2idx) - der_int_T_p1n2w;
            J_TT(rw,elInt_p2n1idx) = J_TT(rw,elInt_p2n1idx) - der_int_T_p2n1w;
            J_TT(rw,elInt_p2n2idx) = J_TT(rw,elInt_p2n2idx) - der_int_T_p2n2w;
            J_TU(rw,elInt_p1n1idx) = J_TU(rw,elInt_p1n1idx) + der_int_T_p1n1w;
            J_TU(rw,elInt_p1n2idx) = J_TU(rw,elInt_p1n2idx) + der_int_T_p1n2w;
            J_TU(rw,elInt_p2n1idx) = J_TU(rw,elInt_p2n1idx) + der_int_T_p2n1w;
            J_TU(rw,elInt_p2n2idx) = J_TU(rw,elInt_p2n2idx) + der_int_T_p2n2w;

            %. third integral: [phi].[u]
            int_A = (gam/dh) * secLen * ( 2*phi_p1*uJ_p1 + phi_p1*uJ_p2 + phi_p2*uJ_p1 + 2*phi_p2*uJ_p2 ) * (1/6);
            F_U(ii,:) = F_U(ii,:) - int_A.';
            F_T(ii,:) = F_T(ii,:) + int_A.';

            %. third integral, deriv.
            der_int_A_p1n1 = (gam/dh) * secLen * ( 2*phi_p1 + phi_p2 ) * elInt_p1n1w * (1/6);
            der_int_A_p1n2 = (gam/dh) * secLen * ( 2*phi_p1 + phi_p2 ) * elInt_p1n2w * (1/6);
            der_int_A_p2n1 = (gam/dh) * secLen * ( phi_p1 + 2*phi_p2 ) * elInt_p2n1w * (1/6);
            der_int_A_p2n2 = (gam/dh) * secLen * ( phi_p1 + 2*phi_p2 ) * elInt_p2n2w * (1/6);
            ind_p1n1 = sub2ind(sz,rw,elInt_p1n1idx);
            ind_p1n2 = sub2ind(sz,rw,elInt_p1n2idx);
            ind_p2n1 = sub2ind(sz,rw,elInt_p2n1idx);
            ind_p2n2 = sub2ind(sz,rw,elInt_p2n2idx);
            J_UU(ind_p1n1(:)) = J_UU(ind_p1n1(:)) + der_int_A_p1n1;
            J_UU(ind_p1n2(:)) = J_UU(ind_p1n2(:)) + der_int_A_p1n2;
            J_UU(ind_p2n1(:)) = J_UU(ind_p2n1(:)) + der_int_A_p2n1;
            J_UU(ind_p2n2(:)) = J_UU(ind_p2n2(:)) + der_int_A_p2n2;
            J_UT(ind_p1n1(:)) = J_UT(ind_p1n1(:)) - der_int_A_p1n1;
            J_UT(ind_p1n2(:)) = J_UT(ind_p1n2(:)) - der_int_A_p1n2;
            J_UT(ind_p2n1(:)) = J_UT(ind_p2n1(:)) - der_int_A_p2n1;
            J_UT(ind_p2n2(:)) = J_UT(ind_p2n2(:)) - der_int_A_p2n2;
            J_TT(ind_p1n1(:)) = J_TT(ind_p1n1(:)) + der_int_A_p1n1;
            J_TT(ind_p1n2(:)) = J_TT(ind_p1n2(:)) + der_int_A_p1n2;
            J_TT(ind_p2n1(:)) = J_TT(ind_p2n1(:)) + der_int_A_p2n1;
            J_TT(ind_p2n2(:)) = J_TT(ind_p2n2(:)) + der_int_A_p2n2;
            J_TU(ind_p1n1(:)) = J_TU(ind_p1n1(:)) - der_int_A_p1n1;
            J_TU(ind_p1n2(:)) = J_TU(ind_p1n2(:)) - der_int_A_p1n2;
            J_TU(ind_p2n1(:)) = J_TU(ind_p2n1(:)) - der_int_A_p2n1;
            J_TU(ind_p2n2(:)) = J_TU(ind_p2n2(:)) - der_int_A_p2n2;
        end
    end

    %. stablisation integrals
    for kk=1:size(inde,2)

        elInd_loc_R = elem_edeges(inde(kk),1);
        elInd_loc_L = elem_edeges(inde(kk),2);
        elInd_R = elem_inds( elInd_loc_R );
        elInd_L = elem_inds( elInd_loc_L );
        elType_R = elTypes( elInd_R );
        elType_L = elTypes( elInd_L );

        if ( elType_R == 3 )||( elType_L == 3 )
            %. one of elements - intersected, stabilise edge

            normInd = elem_edeges(inde(kk),3);

            omeg = 1;
            if ( elem_edeges(inde(kk),1) < 6.5 )&&( elem_edeges(inde(kk),2) < 6.5 )
                %. node ii belongs to the common edge
                omeg = -1;
            end

            nphi = normals(normInd,:).';

            elem_R = elemNodeInds_elem( elInd_R, : );
            elem_L = elemNodeInds_elem( elInd_L, : );

            dGdu_R = elem_R(4) * dGdu_ref * (1/dh);
            dGdu_L = elem_L(4) * dGdu_ref * (1/dh);

            rows = repmat(dofs*ii-sfts.',1,3*dofs);
            idx_R = [ dofs*elem_R(1)-sfts  dofs*elem_R(2)-sfts  dofs*elem_R(3)-sfts ];
            idx_L = [ dofs*elem_L(1)-sfts  dofs*elem_L(2)-sfts  dofs*elem_L(3)-sfts ];
            cols_R = repmat(idx_R,dofs,1);
            cols_L = repmat(idx_L,dofs,1);
            ind_R = sub2ind(sz,rows,cols_R);
            ind_L = sub2ind(sz,rows,cols_L);

            der_int_R = kap * dh * kron(nphi.',eye(dofs)) * dGdu_R * omeg * lengths(normInd)^2;
            der_int_L = kap * dh * kron(nphi.',eye(dofs)) * dGdu_L * omeg * lengths(normInd)^2;

            if ( elType_R ~= 1 )&&( elType_L ~= 1 )
                %. both elements not transformed, stabilise untransformed

                G_U_R = reshape( G_U_elem( elInd_R, : ), dofs, 2 );
                G_U_L = reshape( G_U_elem( elInd_L, : ), dofs, 2 );
                jumpG_U = G_U_R - G_U_L;
                int_U = kap * dh * jumpG_U * nphi * omeg * lengths(normInd)^2;

                F_U(ii,:) = F_U(ii,:) + int_U.';

                J_UU(ind_R(:)) = J_UU(ind_R(:)) + der_int_R(:);
                J_UU(ind_L(:)) = J_UU(ind_L(:)) - der_int_L(:);
            end
            if ( elType_R ~= 2 )&&( elType_L ~= 2 )
                %. both elements not untransformed, stabilise transformed

                G_T_R = reshape( G_T_elem( elInd_R, : ), dofs, 2 );
                G_T_L = reshape( G_T_elem( elInd_L, : ), dofs, 2 );
                jumpG_T = G_T_R - G_T_L;
                int_T = kap * dh * jumpG_T * nphi * omeg * lengths(normInd)^2;

                F_T(ii,:) = F_T(ii,:) + int_T.';

                J_TT(ind_R(:)) = J_TT(ind_R(:)) + der_int_R(:);
                J_TT(ind_L(:)) = J_TT(ind_L(:)) - der_int_L(:);
            end
        end
    end
end

F_U_t = F_U.';
F_T_t = F_T.';
F = [ F_U_t(:); F_T_t(:); ];

J = [ J_UU J_UT ;
      J_TU J_TT ];

%% calculate energy

for ii=1:(2*(Nx-1)*(Ny-1))
    %. volumetric

    W_U = W_U_elem(ii,:);
    W_T = W_T_elem(ii,:);
    if ( elTypes(ii) == 2 )
        %. untransformed
        area_U = 1;
        area_T = 0;
    elseif ( elTypes(ii) == 1 )
        %. transformed
        area_U = 0;
        area_T = 1;
    else
        %. intersected
        indEl_loc = find( intElem == ii, 1, 'first' );
        area_U = fracElem(indEl_loc);
        area_T = 1 - fracElem(indEl_loc);
    end
    E = E + ( area_U * W_U * pweight.' + area_T * W_T * pweight.' ) * (1/2) * (dh^2);
end

for ii=1:size(intElem,1)
    %. interface

    indEl_loc = intElem(ii);
    dWdG_U = reshape( dWdG_U_elem( indEl_loc, : ), dofs, 8 );
    dWdG_T = reshape( dWdG_T_elem( indEl_loc, : ), dofs, 8 );

    dWdG_U_c = dWdG_U(:,1:2);
    dWdG_T_c = dWdG_T(:,1:2);

    elIntPts = intConn(ii,:);
    elInt_p1n1ind = intNodes( elIntPts(1), 1 );
    elInt_p1n1w = fracFaces( elIntPts(1) );
    elInt_p1n2ind = intNodes( elIntPts(1), 2 );
    elInt_p1n2w = 1 - fracFaces( elIntPts(1) );
    elInt_p2n1ind = intNodes( elIntPts(2), 1 );
    elInt_p2n1w = fracFaces( elIntPts(2) );
    elInt_p2n2ind = intNodes( elIntPts(2), 2 );
    elInt_p2n2w = 1 - fracFaces( elIntPts(2) );

    elInt_p1n1idx = dofs*elInt_p1n1ind-sfts;
    elInt_p1n2idx = dofs*elInt_p1n2ind-sfts;
    elInt_p2n1idx = dofs*elInt_p2n1ind-sfts;
    elInt_p2n2idx = dofs*elInt_p2n2ind-sfts;

    uU_p1n1 = u_U( elInt_p1n1idx );
    uU_p1n2 = u_U( elInt_p1n2idx );
    uU_p2n1 = u_U( elInt_p2n1idx );
    uU_p2n2 = u_U( elInt_p2n2idx );
    uU_p1 = uU_p1n1 * elInt_p1n1w + uU_p1n2 * elInt_p1n2w;
    uU_p2 = uU_p2n1 * elInt_p2n1w + uU_p2n2 * elInt_p2n2w;

    uT_p1n1 = u_T( elInt_p1n1idx );
    uT_p1n2 = u_T( elInt_p1n2idx );
    uT_p2n1 = u_T( elInt_p2n1idx );
    uT_p2n2 = u_T( elInt_p2n2idx );
    uT_p1 = uT_p1n1 * elInt_p1n1w + uT_p1n2 * elInt_p1n2w;
    uT_p2 = uT_p2n1 * elInt_p2n1w + uT_p2n2 * elInt_p2n2w;

    uJ_p1 = ( uT_p1 - uU_p1 );
    uJ_p2 = ( uT_p2 - uU_p2 );

    secLen = cutLen(ii);
    secNorm = elemNorms(ii,:).';

    E = E - secLen * ( uJ_p1.' + uJ_p2.' )*(1/2) * ( dWdG_U_c + dWdG_T_c )*(1/2) * secNorm + (gam/dh)*(1/2) * secLen * ( uJ_p1.'*uJ_p1 + uJ_p1.'*uJ_p2 + uJ_p2.'*uJ_p2 )*(1/3);
end

for ii=1:((Nx-1)*(Ny-1))
    %. stablisation

    for jj=1:3
        stabCond = 1;
        if ( jj == 1 )
            %. vertical
            elInd_R = 2*ii + 1;
            elInd_L = 2*ii;
            if ( rem(ii,(Nx-1)) == 0 )
                stabCond = 0;
            end
        elseif ( jj == 2 )
            %. horizontal
            elInd_R = 2*(Nx-1) + 2*ii - 1;
            elInd_L = 2*ii;
            if ( ii > (Nx-1)*(Ny-2) )
                stabCond = 0;
            end
        else
            %. diagonal
            elInd_R = 2*ii;
            elInd_L = 2*ii - 1;
        end

        if ( stabCond == 1 )
            elType_R = elTypes( elInd_R );
            elType_L = elTypes( elInd_L );
            if ( elType_R == 3 )||( elType_L == 3 )
                %. one of elements - intersected, stabilise edge
                nrm = normals(jj,:).';
                len = lengths(jj);
                if ( elType_R ~= 1 )&&( elType_L ~= 1 )
                    %. both elements not transformed, stabilise untransformed
                    G_U_R = reshape( G_U_elem( elInd_R, : ), dofs, 2 );
                    G_U_L = reshape( G_U_elem( elInd_L, : ), dofs, 2 );
                    jumpG_U = G_U_R - G_U_L;
                    E = E + (kap/2) * dh * nrm.' * ( jumpG_U.' * jumpG_U ) * nrm * len * dh;
                end
                if ( elType_R ~= 2 )&&( elType_L ~= 2 )
                    %. both elements not untransformed, stabilise transformed
                    G_T_R = reshape( G_T_elem( elInd_R, : ), dofs, 2 );
                    G_T_L = reshape( G_T_elem( elInd_L, : ), dofs, 2 );
                    jumpG_T = G_T_R - G_T_L;
                    E = E + (kap/2) * dh * nrm.' * ( jumpG_T.' * jumpG_T ) * nrm * len * dh;
                end
            end
        end
    end
end

%% prescribe unused nodes

exclInd_U_r = dofs*repmat(exclNodes_U,1,dofs) - repmat(sfts,size(exclNodes_U,1),1);
exclInd_U = exclInd_U_r(:);
exclInd_T_r = dofs*Nx*Ny + dofs*repmat(exclNodes_T,1,dofs) - repmat(sfts,size(exclNodes_T,1),1);
exclInd_T = exclInd_T_r(:);
F(exclInd_U) = u(exclInd_U);
F(exclInd_T) = u(exclInd_T);
J(exclInd_U,:) = 0;
J(sub2ind(2*sz,exclInd_U,exclInd_U)) = 1;
J(exclInd_T,:) = 0;
J(sub2ind(2*sz,exclInd_T,exclInd_T)) = 1;

%% prescribe coinciding nodes

intIndUnk_U_r = dofs*repmat(intNodesUnk_C,1,dofs) - repmat(sfts,size(intNodesUnk_C,1),1);
intIndUnk_U = intIndUnk_U_r(:);
intIndUnk_T_r = dofs*Nx*Ny + dofs*repmat(intNodesUnk_C,1,dofs) - repmat(sfts,size(intNodesUnk_C,1),1);
intIndUnk_T = intIndUnk_T_r(:);
F(intIndUnk_U) = F(intIndUnk_U) + F(intIndUnk_T);
F(intIndUnk_T) = u(intIndUnk_T) - u(intIndUnk_U);
J(intIndUnk_U,:) = J(intIndUnk_U,:) + J(intIndUnk_T,:);
J(intIndUnk_T,:) = 0;
J(sub2ind(2*sz,intIndUnk_T,intIndUnk_T)) = 1;
J(sub2ind(2*sz,intIndUnk_T,intIndUnk_U)) = -1;

%% boudary conditions

indLB_T = dofs*Nx*Ny+dofs;
indRB_T = dofs*Nx*Ny+dofs*Nx;
indLT_T = dofs*Nx*Ny+dofs*Nx*(Ny-1)+dofs;
indRT_T = dofs*Nx*Ny+dofs*Nx*Ny;

indL_U = (dofs:(dofs*Nx):(dofs*Nx*(Ny-1)+dofs));
indR_U = (dofs:(dofs*Nx):(dofs*Nx*(Ny-1)+dofs))+dofs*(Nx-1);
indB_U = (dofs:dofs:(dofs*Nx));
indT_U = (dofs:dofs:(dofs*Nx))+dofs*Nx*(Ny-1);

indL_T = dofs*Nx*Ny+(dofs:(dofs*Nx):(dofs*Nx*(Ny-1)+dofs));
indR_T = dofs*Nx*Ny+(dofs:(dofs*Nx):(dofs*Nx*(Ny-1)+dofs))+dofs*(Nx-1);
indB_T = dofs*Nx*Ny+(dofs:dofs:(dofs*Nx));
indT_T = dofs*Nx*Ny+(dofs:dofs:(dofs*Nx))+dofs*Nx*(Ny-1);

xx = dh*(0:1:(Nx-1)).';
yy = dh*(0:1:(Ny-1)).';

if ( geomPar.BCtype == 0 )
    %. no rigid body motion, enforced at LB and LT corners (transf.)

    indLBx_T = indLB_T-dofs+1;
    indLBy_T = indLB_T-dofs+2;
    indRBx_T = indRB_T-dofs+1;
    indRBy_T = indRB_T-dofs+2;
    indLTx_T = indLT_T-dofs+1;
    indLTy_T = indLT_T-dofs+2;
    indRTx_T = indRT_T-dofs+1;
    indRTy_T = indRT_T-dofs+2;

    %. corners, displacement
    F(indLBx_T) = u(indLBx_T);
    F(indLBy_T) = u(indLBy_T);
    F(indLTx_T) = u(indLTx_T);

    J( indLBx_T, : ) = 0;
    J( indLBx_T, indLBx_T ) = 1;
    J( indLBy_T, : ) = 0;
    J( indLBy_T, indLBy_T ) = 1;
    J( indLTx_T, : ) = 0;
    J( indLTx_T, indLTx_T ) = 1;

    %. sides, magnetic field
    F(indL_T) = u(indL_T);
    F(indR_T) = u(indR_T);
    F(indB_T) = u(indB_T);
    F(indT_T) = u(indT_T);
    F(indB_U) = u(indB_U);
    F(indT_U) = u(indT_U);

    J( indL_T, : ) = 0;
    J( sub2ind( 2*sz, indL_T, indL_T ) ) = 1;
    J( indR_T, : ) = 0;
    J( sub2ind( 2*sz, indR_T, indR_T ) ) = 1;
    J( indB_T, : ) = 0;
    J( sub2ind( 2*sz, indB_T, indB_T ) ) = 1;
    J( indT_T, : ) = 0;
    J( sub2ind( 2*sz, indT_T, indT_T ) ) = 1;
    J( indB_U, : ) = 0;
    J( sub2ind( 2*sz, indB_U, indB_U ) ) = 1;
    J( indT_U, : ) = 0;
    J( sub2ind( 2*sz, indT_U, indT_U ) ) = 1;

elseif ( geomPar.BCtype == 1 )
    %. L-R uniaxial stretching (transf.)

    displ = geomPar.BCdispl;

    indLx_T = indL_T-dofs+1;
    indLy_T = indL_T-dofs+2;
    indRx_T = indR_T-dofs+1;
    indRy_T = indR_T-dofs+2;
    indBx_T = indB_T-dofs+1;
    indBy_T = indB_T-dofs+2;
    indTx_T = indT_T-dofs+1;
    indTy_T = indT_T-dofs+2;

    indLBy_T = indLB_T-dofs+2;
    
    %. sides, displacement
    F(indLx_T) = u(indLx_T);
    F(indRx_T) = u(indRx_T) - displ;

    J( indLx_T, : ) = 0;
    J( sub2ind( 2*sz, indLx_T, indLx_T ) ) = 1;
    J( indRx_T, : ) = 0;
    J( sub2ind( 2*sz, indRx_T, indRx_T ) ) = 1;

    %. corner, displacement
    F(indLBy_T) = u(indLBy_T);
    
    J( indLBy_T, : ) = 0;
    J( indLBy_T, indLBy_T ) = 1;
    
    %. sides, magnetic field
    F(indL_T) = u(indL_T);
    F(indR_T) = u(indR_T);
    F(indB_T) = u(indB_T);
    F(indT_T) = u(indT_T);
    F(indB_U) = u(indB_U);
    F(indT_U) = u(indT_U);

    J( indL_T, : ) = 0;
    J( sub2ind( 2*sz, indL_T, indL_T ) ) = 1;
    J( indR_T, : ) = 0;
    J( sub2ind( 2*sz, indR_T, indR_T ) ) = 1;
    J( indB_T, : ) = 0;
    J( sub2ind( 2*sz, indB_T, indB_T ) ) = 1;
    J( indT_T, : ) = 0;
    J( sub2ind( 2*sz, indT_T, indT_T ) ) = 1;
    J( indB_U, : ) = 0;
    J( sub2ind( 2*sz, indB_U, indB_U ) ) = 1;
    J( indT_U, : ) = 0;
    J( sub2ind( 2*sz, indT_U, indT_U ) ) = 1;

elseif ( geomPar.BCtype == 2 )
    %. no rigid body motion, enforced at LB and LT corners (transf.)
    %. imposed gradient of magnetic field

    fld = geomPar.BCfld;

    indLBx_T = indLB_T-dofs+1;
    indLBy_T = indLB_T-dofs+2;
    indRBx_T = indRB_T-dofs+1;
    indRBy_T = indRB_T-dofs+2;
    indLTx_T = indLT_T-dofs+1;
    indLTy_T = indLT_T-dofs+2;
    indRTx_T = indRT_T-dofs+1;
    indRTy_T = indRT_T-dofs+2;

    %. corners, displacement
    F(indLBx_T) = u(indLBx_T);
    F(indLBy_T) = u(indLBy_T);
    F(indLTx_T) = u(indLTx_T);

    J( indLBx_T, : ) = 0;
    J( indLBx_T, indLBx_T ) = 1;
    J( indLBy_T, : ) = 0;
    J( indLBy_T, indLBy_T ) = 1;
    J( indLTx_T, : ) = 0;
    J( indLTx_T, indLTx_T ) = 1;

    %. mask for free DOFs
    mskB_T = ~ismember(indB_T,exclInd_T);
    mskT_T = ~ismember(indT_T,exclInd_T);
    mskB_U = ~ismember(indB_U,exclInd_U);
    mskT_U = ~ismember(indT_U,exclInd_U);

    %. L-R sides, magnetic field
    % F(indL_T) = u(indL_T) - fld;
    % F(indR_T) = u(indR_T);

    % J( indL_T, : ) = 0;
    % J( sub2ind( 2*sz, indL_T, indL_T ) ) = 1;
    % J( indR_T, : ) = 0;
    % J( sub2ind( 2*sz, indR_T, indR_T ) ) = 1;

    %. B-T sides, magnetic field
    F(indB_T(mskB_T)) = u(indB_T(mskB_T)) - fld;
    F(indT_T(mskT_T)) = u(indT_T(mskT_T));
    F(indB_U(mskB_U)) = u(indB_U(mskB_U)) - fld;
    F(indT_U(mskT_U)) = u(indT_U(mskT_U));

    J( indB_T(mskB_T), : ) = 0;
    J( sub2ind( 2*sz, indB_T(mskB_T), indB_T(mskB_T) ) ) = 1;
    J( indT_T(mskT_T), : ) = 0;
    J( sub2ind( 2*sz, indT_T(mskT_T), indT_T(mskT_T) ) ) = 1;
    J( indB_U(mskB_U), : ) = 0;
    J( sub2ind( 2*sz, indB_U(mskB_U), indB_U(mskB_U) ) ) = 1;
    J( indT_U(mskT_U), : ) = 0;
    J( sub2ind( 2*sz, indT_U(mskT_U), indT_U(mskT_U) ) ) = 1;

elseif ( geomPar.BCtype == 3 )
    %. sliding L/R/B/T, biaxial stretching (transf.)

    displ = geomPar.BCdispl;

    indLx_T = indL_T-dofs+1;
    indLy_T = indL_T-dofs+2;
    indRx_T = indR_T-dofs+1;
    indRy_T = indR_T-dofs+2;
    indBx_T = indB_T-dofs+1;
    indBy_T = indB_T-dofs+2;
    indTx_T = indT_T-dofs+1;
    indTy_T = indT_T-dofs+2;

    %. B.C. transf. bottom
    F(indBy_T) = u(indBy_T);

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T);

    %. B.C. transf. right
    F(indRx_T) = u(indRx_T) - displ;

    %. B.C. transf. top
    F(indTy_T) = u(indTy_T) - displ;

    J( indBy_T, : ) = 0;
    J( sub2ind( 2*sz, indBy_T, indBy_T ) ) = 1;
    J( indLx_T, : ) = 0;
    J( sub2ind( 2*sz, indLx_T, indLx_T ) ) = 1;
    J( indRx_T, : ) = 0;
    J( sub2ind( 2*sz, indRx_T, indRx_T ) ) = 1;
    J( indTy_T, : ) = 0;
    J( sub2ind( 2*sz, indTy_T, indTy_T ) ) = 1;

elseif ( geomPar.BCtype == 4 )
    %. magnetic domain wall

    thetL = asin(tanh(-pi*sqrt(2)/2)) + pi/2;
    thetR = asin(tanh(pi*sqrt(2)/2)) + pi/2;

    indLx_T = indL_T-dofs+1;
    indLy_T = indL_T-dofs+2;
    indLz_T = indL_T-dofs+3;
    indLl_T = indL_T-dofs+4;
    indRx_T = indR_T-dofs+1;
    indRy_T = indR_T-dofs+2;
    indRz_T = indR_T-dofs+3;
    indRl_T = indR_T-dofs+4;

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T) - sin(thetL);
    F(indLy_T) = u(indLy_T) - cos(thetL);
    F(indLz_T) = u(indLz_T);
    F(indLl_T) = u(indLl_T) - 1;

    %. B.C. transf. right
    F(indRx_T) = u(indRx_T) - sin(thetR);
    F(indRy_T) = u(indRy_T) - cos(thetR);
    F(indRz_T) = u(indRz_T);
    F(indRl_T) = u(indRl_T) - 1;

    J( indLx_T, : ) = 0;
    J( sub2ind( 2*sz, indLx_T, indLx_T ) ) = 1;
    J( indLy_T, : ) = 0;
    J( sub2ind( 2*sz, indLy_T, indLy_T ) ) = 1;
    J( indLz_T, : ) = 0;
    J( sub2ind( 2*sz, indLz_T, indLz_T ) ) = 1;
    J( indLl_T, : ) = 0;
    J( sub2ind( 2*sz, indLl_T, indLl_T ) ) = 1;
    J( indRx_T, : ) = 0;
    J( sub2ind( 2*sz, indRx_T, indRx_T ) ) = 1;
    J( indRy_T, : ) = 0;
    J( sub2ind( 2*sz, indRy_T, indRy_T ) ) = 1;
    J( indRz_T, : ) = 0;
    J( sub2ind( 2*sz, indRz_T, indRz_T ) ) = 1;
    J( indRl_T, : ) = 0;
    J( sub2ind( 2*sz, indRl_T, indRl_T ) ) = 1;

end

end


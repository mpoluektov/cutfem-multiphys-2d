function [ ] = mainPhaseTran( filename )
%MAINPHASETRAN Cut finite elements for modelling phase transitions
%   See (Poluektov and Figiel, Comput. Mech. 63:885-911, 2019).
%   Code has been written by M. Poluektov.

%% param.

if ( nargin < 1 )
    error( 'Specify filename' );
end

%. elastic properties
matPar.elGU = 24;
matPar.elGT = 10;
matPar.elKU = 75;
matPar.elKT = 32;
matPar.trfU = 1;
matPar.trfT = 1.05;
% matPar.elGU = 2;
% matPar.elGT = 2;
% matPar.elKU = 100;
% matPar.elKT = 100;
% matPar.trfU = 1;
% matPar.trfT = 1;

%. magnetic properties
matPar.mgAU = 1;
matPar.mgAT = 1;
matPar.mgKU = 2*pi^2;
matPar.mgKT = 2*pi^2;
matPar.aniU = pi/2;
matPar.aniT = pi/2;
% matPar.mgAU = 1/6000;
% matPar.mgAT = 1/6000;
% matPar.mgKU = 2*pi^2/60;
% matPar.mgKT = 2*pi^2/60;
% matPar.aniU = pi/2;
% matPar.aniT = 0;
matPar.cf1U = 1/100;
matPar.cf1T = 1/100;
matPar.cf2U = 1/1000;
matPar.cf2T = 1/1000;
matPar.cf3U = 0;
matPar.cf3T = 0;

%. chemical properties
matPar.gam = -0.15;
matPar.coefVel = 0.043;
% matPar.gam = 0;
% matPar.coefVel = 0.05;

%. grid
domainSize = 1;
geomPar.Nx = 24 + 1;
geomPar.Ny = 24 + 1;
% geomPar.Nx = 63 + 1;
% geomPar.Ny = 21 + 1;
% geomPar.Nx = 64 + 1;
% geomPar.Ny = 64 + 1;
geomPar.dh = domainSize/(geomPar.Nx-1);

%. number of degrees of freedom
geomPar.dofs = 2;
% geomPar.dofs = 4;

%. boundary conditions for mechanics
% geomPar.BCtype = 0;
% geomPar.BCtype = 1;
% geomPar.BCdispl = 0.02;
% geomPar.BCdispl = -0.07;
% geomPar.BCtype = 2;
% geomPar.BCfld = 4;
% geomPar.BCfld = 4/3;
geomPar.BCtype = 3;
geomPar.BCdispl = 0.038;
% geomPar.BCtype = 4;

%. numerical parameters
numPar.ATOL = 1e-11;
numPar.DTOL = 1e-2;
numPar.MAX_ITER = 30;
numPar.SDLT = 1e-8;
numPar.NIT_PAR = 1e4;
% numPar.NIT_PAR = 1;
numPar.STAB_PAR = 1e-1;
% numPar.STAB_PAR = 1;
numPar.PTOL = 1e-14;

%. time
tend = 100;
Nt = 100 + 1;
% tend = 200;
% Nt = 800 + 1;
% Nt = 1;
dt = tend/(Nt-1);

%% initialise

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;
dofs = geomPar.dofs;

%. create grid
gXs = ( 0:1:(Nx-1) )*dh;
gYs = ( 0:1:(Ny-1) )*dh;
gX = repmat(gXs,1,Ny);
gYa = repmat(gYs,1,Nx);
gYb = reshape(gYa,Ny,[]);
gYt = gYb.';
gY = gYt(:).';
coordGrid = [ gX; gY; ];

%% create initial interface

intType = 1;

if ( intType == 0 )
    %. flat interface
    
    geomPar.windThresh = 0;
    
    y0 = 0.1;
    
    intP = [  -dh y0 ;
             1+dh y0 ];

    conn = [ 1 2 ];

elseif ( intType == 1 )
    %. circular interface

    geomPar.windThresh = 0.5;
    
    x0 = 0.5;
    y0 = 0.5;
    a = 0.27;
    % a = 2/7;

    dalp = pi/180;
    st = 0:360;
    px = x0 + a*cos(dalp*st);
    py = y0 + a*sin(dalp*st);
    intP = [ px.' py.' ];
    
    NintP = size( intP, 1 );
    conn = [ 1:(NintP-1) ;
             2:NintP     ].';

elseif ( intType == 2 )
    %. twin boundary

    geomPar.windThresh = 0.1;

    intP = [ 1/2 0   ;
             1/6 1/3 ;
             5/6 0   ;
             1/2 1/3 ];

    conn = [ 3 4 ;
             2 1 ];

elseif ( intType == 4 )
    %. several circular interfaces

    geomPar.windThresh = 0.5;

    x0 = 0.65;
    y0 = 0.65;
    a = 0.11;

    dalp = pi/180;
    st = 0:360;
    NintP = size( st, 2 );

    px = x0 + a*cos(dalp*st);
    py = y0 + a*sin(dalp*st);
    intP1 = [ px.' py.' ];
    conn1 = [ 1:(NintP-1) ;
              2:NintP     ].';

    x0 = 0.65+0.35*cos(pi/12+pi/3+pi/2);
    y0 = 0.65-0.35*sin(pi/12+pi/3+pi/2);
    a = 0.11;

    px = x0 + a*cos(dalp*st);
    py = y0 + a*sin(dalp*st);
    intP2 = [ px.' py.' ];
    conn2 = [ (NintP+1):(2*NintP-1) ;
              (NintP+2):(2*NintP)   ].';

    x0 = 0.65+0.35*cos(pi/12+pi/2);
    y0 = 0.65-0.35*sin(pi/12+pi/2);
    a = 0.11;

    px = x0 + a*cos(dalp*st);
    py = y0 + a*sin(dalp*st);
    intP3 = [ px.' py.' ];
    conn3 = [ (2*NintP+1):(3*NintP-1) ;
              (2*NintP+2):(3*NintP)   ].';

    intP = [ intP1; intP2; intP3 ];
    conn = [ conn1; conn2; conn3 ];

end

%% plot

plotInt = 0;
if ( plotInt == 1 )
    tri = zeros( 2*(Nx-1)*(Ny-1), 3 );
    for ii=1:((Nx-1)*(Ny-1))

        ind_bl = ii + floor((ii-0.5)/(Nx-1));
        ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
        ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
        ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));

        tri( 2*ii-1, 1 ) = ind_br;
        tri( 2*ii-1, 2 ) = ind_bl;
        tri( 2*ii-1, 3 ) = ind_tl;
        tri( 2*ii, 1 ) = ind_tl;
        tri( 2*ii, 2 ) = ind_tr;
        tri( 2*ii, 3 ) = ind_br;
    end

    figure(1);
    hold on;

    triplot( tri, coordGrid(1,:).', coordGrid(2,:).', 'Color', [0.75 0.75 0.75] );

    plot( intP(:,1), intP(:,2), 'd-' );

    limits = [ -0.1 1.1 -0.1 1.1 ];
    axis( limits );
    pbaspect( [ 1 1 1 ] );
end

%% solve and move front

ppStart = 1;
mkdir( filename );
save( [ filename '/param' ], 'coordGrid', 'matPar', 'geomPar', 'numPar', 'Nt', 'tend', 'dt' );

% ppStart = 80;
% load( [ filename sprintf('/res%05i',ppStart-1) ], 'res' );
% intP = res.intP_new;
% conn = res.intParam.intConn;
% ucmr = res.ucmr;

% load( 'approxFreeInt', 'intP', 'conn' );

for pp=ppStart:Nt

    %. intersection points and elements
    intParam = calcIntLevelSet( intP, conn, geomPar, numPar.PTOL );

    if ( plotInt == 1 )
        %. transf. elem.
        triplot( tri(intParam.elTypes==1,:), coordGrid(1,:).', coordGrid(2,:).', 'Color', 'b' );
        %. untransf. elem.
        triplot( tri(intParam.elTypes==2,:), coordGrid(1,:).', coordGrid(2,:).', 'Color', 'm' );
    end
    
    %. initial estimate
    if ( pp == 1 )
        u0 = zeros( 2*dofs*Nx*Ny, 1 );
        % u0 = rand( 2*dofs*Nx*Ny, 1 )/100;

        % load( 'approxFreeSol', 'u' );
        % u0 = u;
    else
        %. better initial estimate
        urUt = ucmr;
        urTt = ucmr;
        urUt(:,intParam.exclNodes_U) = 0;
        urTt(:,intParam.exclNodes_T) = 0;
        ur = [ urUt urTt ];
        u0 = ur(:);
    end
    
    %. initial estimate for magnetic domain wall
    % ucmr = zeros(4,Nx*Ny);
    % ucmr(1,:) = sin( 0.2 + (pi-0.4)*coordGrid(1,:) );
    % ucmr(2,:) = cos( 0.2 + (pi-0.4)*coordGrid(1,:) );
    % ucmr(4,:) = 1;
    % urUt = ucmr;
    % urTt = ucmr;
    % urUt(:,intParam.exclNodes_U) = 0;
    % urTt(:,intParam.exclNodes_T) = 0;
    % ur = [ urUt urTt ];
    % u0 = ur(:);

    %. solve one timestep
    [ u, noConv, resElem ] = solveNR( u0, matPar, geomPar, numPar, intParam );
    if ( noConv ~= 0 )
        break;
    end
    
    %. reshape solution
    ur = reshape(u,dofs,2*Nx*Ny);
    urUt = ur(:,1:(Nx*Ny));
    urTt = ur(:,(Nx*Ny+1):(2*Nx*Ny));
    urUt(:,intParam.exclNodes_U) = 0;
    urUt(:,intParam.intNodesUnk_U) = 0;
    urUt(:,intParam.intNodesUnk_C) = urUt(:,intParam.intNodesUnk_C)/2;
    urTt(:,intParam.exclNodes_T) = 0;
    urTt(:,intParam.intNodesUnk_T) = 0;
    urTt(:,intParam.intNodesUnk_C) = urTt(:,intParam.intNodesUnk_C)/2;
    ucmr = urUt + urTt;
    
    %. solution at the interface
    NintP = size(intParam.intPoints,1);
    ui = zeros(NintP,dofs);
    for ii=1:NintP
        n1 = intParam.intNodes(ii,1);
        n2 = intParam.intNodes(ii,2);
        w1 = intParam.fracFaces(ii);
        w2 = 1 - intParam.fracFaces(ii);
        un1 = ucmr(:,n1).';
        un2 = ucmr(:,n2).';
        ui(ii,:) = un1*w1 + un2*w2 ;
    end
    
    %. energy at the interface
    [ resAv ] = interElemAver( resElem, geomPar, intParam );
    [ WUi, WTi, GUi, GTi, jmpi ] = intEnerg( resAv, geomPar, intParam );
    xi = matPar.gam - WTi + WUi + jmpi;

    %. kinetics without field for magneto-mechanics
    % xi_mm = zeros(NintP,1);
    % for ii=1:NintP
    %     [ W_U, dWdG_U ] = constLawCoupledPlan( GUi(ii,:), ui(ii,:), matPar.elKU, matPar.elGU, matPar.trfU, matPar.mgAU, matPar.mgKU, matPar.aniU, 0, 0, 0, numPar.SDLT );
    %     [ W_T, dWdG_T ] = constLawCoupledPlan( GTi(ii,:), ui(ii,:), matPar.elKT, matPar.elGT, matPar.trfT, matPar.mgAT, matPar.mgKT, matPar.aniT, 0, 0, 0, numPar.SDLT );
    %     jmp = ( dWdG_T + dWdG_U ) * (1/2) * ( GTi(ii,:) - GUi(ii,:) ).';
    %     xi_mm(ii) = matPar.gam - W_T + W_U + jmp;
    % end

    %. move interface
    intP_new = moveInt( xi, matPar, intParam, dt );
    % intP_new = moveInt( xi_mm, matPar, intParam, dt );

    if ( plotInt == 1 )
        plot( intP_new(:,1), intP_new(:,2), 'd' );
    end

    %. export
    res.u = u;
    res.ucmr = ucmr;
    res.u_int = ui;
    res.resElem = resElem;
    res.resAv = resAv;
    res.intParam = intParam;
    res.WUi = WUi;
    res.WTi = WTi;
    res.jmpi = jmpi;
    res.intP_new = intP_new;

    %. update interface points
    intP = intP_new;
    conn = intParam.intConn;
    
    fprintf( '  step %.0f/%.0f done\n', pp, Nt );
    
    %. save
    save( [ filename '/res' sprintf('%05i',pp) ], 'res' );

end

end


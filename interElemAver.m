function [ resAv ] = interElemAver( resElem, geomPar, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dofs = geomPar.dofs;

W_U_elem = resElem.W_U_elem(:,1);
W_T_elem = resElem.W_T_elem(:,1);
G_U_elem = resElem.G_U_elem;
G_T_elem = resElem.G_T_elem;
dWdG_U_elem = resElem.dWdG_U_elem(:,1:(2*dofs));
dWdG_T_elem = resElem.dWdG_T_elem(:,1:(2*dofs));

W_U_av = zeros( Nx*Ny, 1 );
W_T_av = zeros( Nx*Ny, 1 );
G_U_av = zeros( Nx*Ny, 2*dofs );
G_T_av = zeros( Nx*Ny, 2*dofs );
dWdG_U_av = zeros( Nx*Ny, 2*dofs );
dWdG_T_av = zeros( Nx*Ny, 2*dofs );

intElem = intParam.intElem;
fracElem = intParam.fracElem;
elTypes = intParam.elTypes;

for ii=1:(Nx*Ny)
               
    %. squares
    sq_ar = ii - floor( (ii-0.5)/Nx );
    sq_al = ii - floor( (ii-0.5)/Nx ) - 1;
    sq_br = ii - floor( (ii-0.5)/Nx ) - (Nx-1);
    sq_bl = ii - floor( (ii-0.5)/Nx ) - (Nx-1) - 1;

    elem_inds = [ 2*sq_al-1 ;
                  2*sq_ar-1 ;
                  2*sq_br-1 ;
                    2*sq_al ;
                    2*sq_br ;
                    2*sq_bl ];

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

    area_U = zeros( size(inda,2), 1 );
    area_T = zeros( size(inda,2), 1 );
    W_U = zeros( size(inda,2), 1 );
    W_T = zeros( size(inda,2), 1 );
    G_U = zeros( size(inda,2), 2*dofs );
    G_T = zeros( size(inda,2), 2*dofs );
    dWdG_U = zeros( size(inda,2), 2*dofs );
    dWdG_T = zeros( size(inda,2), 2*dofs );
    
    for jj=1:size(inda,2)
        
        elType_loc = elTypes( elem_inds(inda(jj)) );
        if ( elType_loc == 2 )
            %. untransformed
            area_U(jj) = 1;
            area_T(jj) = 0;
        elseif ( elType_loc == 1 )
            %. transformed
            area_U(jj) = 0;
            area_T(jj) = 1;
        else
            %. intersected
            useFracArea = 1;
            if ( useFracArea == 1 )
                indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
                area_U(jj) = fracElem(indEl_loc);
                area_T(jj) = 1 - fracElem(indEl_loc);
            else
                area_U(jj) = 1;
                area_T(jj) = 1;
            end
        end
        
        W_U(jj,:) = W_U_elem( elem_inds(inda(jj)), : );
        W_T(jj,:) = W_T_elem( elem_inds(inda(jj)), : );
        G_U(jj,:) = G_U_elem( elem_inds(inda(jj)), : );
        G_T(jj,:) = G_T_elem( elem_inds(inda(jj)), : );
        dWdG_U(jj,:) = dWdG_U_elem( elem_inds(inda(jj)), : );
        dWdG_T(jj,:) = dWdG_T_elem( elem_inds(inda(jj)), : );
    end
    
    area_U_sum = sum(area_U);
    area_T_sum = sum(area_T);
    if ( area_U_sum ~= 0 )
        W_U_av(ii,:) = area_U.' * W_U / area_U_sum;
        G_U_av(ii,:) = area_U.' * G_U / area_U_sum;
        dWdG_U_av(ii,:) = area_U.' * dWdG_U / area_U_sum;
    end
    if ( area_T_sum ~= 0 )
        W_T_av(ii,:) = area_T.' * W_T / area_T_sum;
        G_T_av(ii,:) = area_T.' * G_T / area_T_sum;
        dWdG_T_av(ii,:) = area_T.' * dWdG_T / area_T_sum;
    end
end

resAv.W_U_av = W_U_av;
resAv.W_T_av = W_T_av;
resAv.G_U_av = G_U_av;
resAv.G_T_av = G_T_av;
resAv.dWdG_U_av = dWdG_U_av;
resAv.dWdG_T_av = dWdG_T_av;

end


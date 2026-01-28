function [ u, noConv, resElem ] = solveNR( u0, matPar, geomPar, numPar, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dofs = geomPar.dofs;

ATOL = numPar.ATOL;
DTOL = numPar.DTOL;
MAX_ITER = numPar.MAX_ITER;
SDLT = numPar.SDLT;

u = u0;

noConv = 1;
nrIter = 1;
while ( noConv == 1 )
    
    [ E, F, J, resElem ] = calcRes( u, matPar, geomPar, numPar, intParam );
    
    %. numerical Jacobian
    calcNumJ = 0;
    if ( calcNumJ == 1 )
        F_num = zeros(2*dofs*Nx*Ny,1);
        J_num = zeros(2*dofs*Nx*Ny,2*dofs*Nx*Ny);
        for kk = 1:(2*dofs*Nx*Ny)
            u_c = u;
            u_c(kk) = u_c(kk) + SDLT;
            [ E_c, F_c ] = calcRes( u_c, matPar, geomPar, numPar, intParam );
            F_num( kk ) = ( E_c - E ) / SDLT;
            J_num( :, kk ) = ( F_c - F ) / SDLT;
            disp(kk);
        end
    end

    F_s = sparse(F);
    J_s = sparse(J);
    Dlt_s = -J_s \ F_s;
    Dlt = full(Dlt_s);

    normF = max(abs(F));
    normDlt = max(abs(Dlt));
    
    regSol = 0;
    lamb = 1;
    if ( max(normF,normDlt) > DTOL )&&( regSol == 1 )
        lamb = 0.5;
    end

    u = u + lamb*Dlt;
    
    if ( normF < ATOL )&&( normDlt < ATOL )&&( isreal(Dlt) )
        noConv = 0;
        fprintf( '    N.R. converged, %.0f iter., normF=%1.1e, normDlt=%1.1e\n', nrIter, normF, normDlt );
    else
        fprintf( '        iteration completed, normF=%1.1e, normDlt=%1.1e, lamb=%0.1f, isreal(Dlt)=%0.0f\n', normF, normDlt, lamb, isreal(Dlt) );
    end
    if ( nrIter > MAX_ITER )
        noConv = 2;
        fprintf( '    NO CONVERGENCE, N.R.\n' );
    end
    if ( ~isreal(Dlt) )
        noConv = 3;
        fprintf( '    NO CONVERGENCE, N.R.\n' );
    end
    nrIter = nrIter + 1;
    
end

end


function [ intPoints_new ] = moveInt( xi, matPar, intParam, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

coefVel = matPar.coefVel;

intPoints = intParam.intPoints;
intNorms = intParam.intNorms;

%. velocities, phase transition
vel_int = coefVel * xi;

%. new positions of interface points
intPoints_new = intPoints + dt * intNorms .* repmat( vel_int, 1, 2 );

end


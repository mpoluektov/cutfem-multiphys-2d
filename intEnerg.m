function [ WUi, WTi, GUi, GTi, jmpi ] = intEnerg( resAv, geomPar, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dofs = geomPar.dofs;

fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;

NintP = size(intNodes,1);

WUi = zeros(NintP,1);
WTi = zeros(NintP,1);
GUi = zeros(NintP,2*dofs);
GTi = zeros(NintP,2*dofs);
jmpi = zeros(NintP,1);

W_U_av = resAv.W_U_av;
W_T_av = resAv.W_T_av;
G_U_av = resAv.G_U_av;
G_T_av = resAv.G_T_av;
dWdG_U_av = resAv.dWdG_U_av;
dWdG_T_av = resAv.dWdG_T_av;

for ii=1:NintP

    n1 = intNodes(ii,1);
    n2 = intNodes(ii,2);
    w1 = fracFaces(ii);
    w2 = 1 - fracFaces(ii);

    WU1 = W_U_av(n1,:);
    WT1 = W_T_av(n1,:);
    WU2 = W_U_av(n2,:);
    WT2 = W_T_av(n2,:);
    GU1 = G_U_av(n1,:);
    GT1 = G_T_av(n1,:);
    GU2 = G_U_av(n2,:);
    GT2 = G_T_av(n2,:);
    dWdGU1 = dWdG_U_av(n1,:);
    dWdGT1 = dWdG_T_av(n1,:);
    dWdGU2 = dWdG_U_av(n2,:);
    dWdGT2 = dWdG_T_av(n2,:);
    
    WUa = WU1*w1 + WU2*w2;
    WTa = WT1*w1 + WT2*w2;
    GUa = GU1*w1 + GU2*w2;
    GTa = GT1*w1 + GT2*w2;
    dWdGUa = dWdGU1*w1 + dWdGU2*w2;
    dWdGTa = dWdGT1*w1 + dWdGT2*w2;
    
    jmp = ( dWdGTa + dWdGUa ) * (1/2) * ( GTa - GUa ).';
    
    WUi(ii) = WUa;
    WTi(ii) = WTa;
    GUi(ii,:) = GUa;
    GTi(ii,:) = GTa;
    jmpi(ii) = jmp;

end

end


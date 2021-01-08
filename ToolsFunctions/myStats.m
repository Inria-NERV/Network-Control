%{
this function computes all predefined stats in the following order:
--------- Trajectory -------
success
energy cost
trajectory length
-------- Controllability ----
kalman rank
grammian rank
grammian condition number
average controllability
PBH test

%}

function [success,en,L,rC , rW, condGram, TrGram, pbhTest]= myStats(x , u,x0,xf, eps , A,B,t,STEP )
%% Trajectory
success = controlSucces(x,xf,x0, eps);
en = controlEnergy(u);
L = trajLength(x);
%% Controllability
C = ctrb(A,B);
rC = rank(C);
W = myGram(A,B,t,STEP);
rW = rank(W);
TrGram = trace(W);
condGram = log10(cond(W));
pbhTest = myPBHtest(A,B);
%% summary
%stats = [ s en L rC rW condGram TrGram pbhTest] ;
end
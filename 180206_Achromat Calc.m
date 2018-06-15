clear all
%%
% Initial Guess of Glass: N-BK7 (Crown, convex), N-SF10 (Flint, concave)
F_Target=90;    %mm
% Choose N s-line (852nm) and V d-line (587.6nm)
N1=1.50980;%s-line: 1.50980;
V1=65;%d-line: 64.17;
N2=1.70891;%s-line: 1.70891;
V2=29;%d-line: 28.53;

R1=F_Target*(V1-V2)/V1*(N1-1)

R2=F_Target*(V1-V2)/V2*(N2-1)


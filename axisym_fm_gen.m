clear
clc

syms lnkf lnCf lnaf lnhf lnks_perp lnks_par lnCs lnas lnRth lnsx lnsy f u v vs1 vs2 vs3

ks_perp = exp(lnks_perp);
Ks = ks_perp * eye(3) + (exp(lnks_par) - ks_perp) * [vs1;vs2;vs3] * [vs1,vs2,vs3];

Kf11 = exp(lnkf);
Kf12 = 0;
Kf13 = 0;
Kf22 = exp(lnkf);
Kf23 = 0;
Kf33 = exp(lnkf);
Cf   = exp(lnCf);
af   = exp(lnaf);
Rf   = 0;
hf   = exp(lnhf);
Ks11 = Ks(1,1);
Ks12 = Ks(1,2);
Ks13 = Ks(1,3);
Ks22 = Ks(2,2);
Ks23 = Ks(2,3);
Ks33 = Ks(3,3);
Cs   = exp(lnCs);
as   = exp(lnas);
Rs   = 0;
hs   = 1000;
Rth  = exp(lnRth);
sx   = exp(lnsx);
sy   = exp(lnsy);
P    = 1;
f0   = f;

T0hat = T0hat_infinite(Kf11,Kf12,Kf13,Kf22,Kf23,Kf33,Cf,af,Rf,hf,Ks11,Ks12,Ks13,Ks22,Ks23,Ks33,Cs,as,Rs,hs,Rth,sx,sy,P,f0,u,v);
grad = gradient(T0hat, [lnkf lnCf lnhf lnks_perp lnks_par lnCs vs1 vs2 vs3]);

m_str = ["lnkf" "lnCf" "lnhf" "lnks_perp" "lnks_par" "lnCs" "v1" "v2" "v3"];
g_str = "g_" + m_str;

in_vars = [lnkf lnCf lnaf lnhf lnks_perp lnks_par lnCs lnas lnRth lnsx lnsy f u v vs1 vs2 vs3];

matlabFunction(T0hat, Optimize=true, File="axisym_fm.m", Vars=in_vars, Outputs="T0hat")
matlabFunction(T0hat, grad(1), grad(2), grad(3), grad(4), grad(5), grad(6), grad(7), grad(8), grad(9), ...
    Optimize=true, File="axisym_fm_g.m", Vars=in_vars, ...
    Outputs=["T0hat", g_str(:).'])
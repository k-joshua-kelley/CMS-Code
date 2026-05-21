clear
clc

syms lnkf lnCf lnaf lnhf lnks lnCs lnas lnRth lnsx lnsy f u v

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
Ks11 = exp(lnks);
Ks12 = 0;
Ks13 = 0;
Ks22 = exp(lnks);
Ks23 = 0;
Ks33 = exp(lnks);
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
grad = gradient(T0hat, [lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth]);

m_str = ["lnkf", "lnCf", "lnaf", "lnhf", "lnks", "lnCs", "lnas", "lnRth"];
g_str = "g_" + m_str;

in_vars = [lnkf lnCf lnaf lnhf lnks lnCs lnas lnRth lnsx lnsy f u v];

matlabFunction(T0hat, Optimize=true, File="iso_fm.m", Vars=in_vars, Outputs="T0hat")
matlabFunction(T0hat, grad(1), grad(2), grad(3), grad(4), grad(5), grad(6), grad(7), grad(8), ...
    Optimize=true, File="iso_fm_g.m", Vars=in_vars, ...
    Outputs=["T0hat", g_str(:).'])
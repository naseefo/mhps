from mhps.isolator_osbi import fb, fbfixed, fs1, fs1fixed, fs2, fs2fixed, IsoOSBIModel

rmbm = 1.0
tbx = 50.0
zetabx = 0.0
rtytxb = 1.0
rzyzxb = 1.0
typevf = 0
mu0 = 0.1
alpha0 = 1
alpha1 = 1
nu = 1
umax = 0.85
D = 0.5
ecc = 0.4
rmrm = 0.05
tc = 100.0
ecrc = 0.0
fos_ud = 1.0
am = 1.0
niso = 4.0

iso = IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, tc, ecrc, fos_ud, am, niso)
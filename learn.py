# import click
# import re

# @click.command()
# @click.option('--eq','-e', default='db',type=str)
# def main(eq):
#     eq_1 = re.findall("Excitations.csv",eq)
#     eq_2 = re.findall("\[[\w\-\.]+.[txt|csv], [\w\-\.]+.[txt|csv], [0-9\.]+, [\w/]+, [0-9\.]+, [a-z]]",eq)
#     eq_3 = re.findall("\[[\w\-\.]+.[txt|csv], [0-9\.]+, [\w/]+, [0-9\.]+, [a-z]]",eq)
#     eq_4 = re.findall("db",eq)

#     print(len(eq_1))
#     print(len(eq_2))
#     print(len(eq_3))
#     print(len(eq_4))

#     for item in eq_1:
#         print(item)

#     for item in eq_2:
#         print(item)

#     for item in eq_3:
#         print(item)

#     for item in eq_4:
#         print(item)
    



# if __name__ == '__main__':
#     main()


# # Excitations.csv
# [Cent_acc_00, Cent_acc_00, 0.02, cm/s2, 20]
# [eq-x1.txt, eq.wwy1.txt, 0.02, cm/s2, 20], [eqx2.txt, eqy2.txt, 0.005, g, 30]
# [eqx1.rt, eq.wwy1.txt, 0.02, cm/s2, 20], [eqx2.txt, eqy2.txt, 0.005, g, 30]

# [eqx1.txt, 0.02, cm/s2, 20, x], [eqx2.txt, 0.005, g, 30, y]

# # [eqx1.txt], [eqx2.txt]

# # [eqx1.txt eqx1.txt], [eqx2.txt eqx1.txt]

# Cent_acc_00.txt, Cent_acc_00.txt, 20;


# Cent_acc_00.txt, 10, x


[AXI, TX1, PHX1, AY1, TY1, PHY1, DT, UNIT, DUR, NDIV]
[3.0, 0.1, 1.5, 2.0, 0.2, 1.0, 0.02, cm/s2, 20]
\[[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[\w/]+[,\s]+[0-9\.]+]

[AXI, TX1, PHX1, DT, UNIT, DUR, DIR]
[3.0, 0.1, 1.5, 0.02, cm/s2, 20, x]
[3.0, 0.1, 1.5, 0.02, cm/s2, 20, y]
\[[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[\w/]+[,\s]+[0-9\.]+[,\s]+[\w]]

[AX1, TX1_0:DTX1:TX1_N, PHX1, DT, UNIT, DUR, DIR]
[3.0, 0.1:0.01:5.0, 1.5, 0.02, cm/s2, 20, x]
\[[0-9\.]+[,\s]+[0-9\.]+:[0-9\.]+:[0-9\.]+[,\s]+[0-9\.]+[,\s]+[0-9\.]+[,\s]+[\w/]+[,\s]+[0-9\.]+[,\s]+[\w]]




# import glob, os

# # os.chdir("/")


# files = glob.glob('*.txt')
# files.extend(glob.glob('Excitations.csv'))
# print(files)


# def store_sq(num):
#     for i in range(0,num):
#         yield i^2

^([\w\.\s]+),([\w\.\s]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9]+),([\s0-9]+)(,[\s0-9]+|)\;

Cent_acc_00.txt, Cent_acc_90.txt, 0.02, cm/s2, 1, 30, 100;

^([\w\.\s]+)[,\s]+([s0-9\.]+)[,\s]+([a-z0-9/]+)[,\s]([\s0-9\.]+)[,\s]+([0-9]+)[,\s]+([\w])[,\s]+([0-9]+);
Cent_acc_00.csv, 0.02, g, 1.0, 10, y, 120;
Cent_acc_00.csv, 0.02, g, 10, y, 120;



* - All gx, dx, vx, ax, aax, gy, dy, vy, ay, aay, fx, fy, ek, ed, es, ei, error



      - gx, dx, vx, ax, aax, gy, dy, vy, ay, aay, fx, fy, ek, ed, es, ei, error
      - gx, gy, dx, vx, ax, aax, dy, vy, ay, aay



aa*
a*
v*
d*
r*

aa2
a2
v2
d2
r2

aa3-5
a3-5
v3-5
d3-5
r3-5

aab
ab
vb
db
rb

paa2

f
e
er


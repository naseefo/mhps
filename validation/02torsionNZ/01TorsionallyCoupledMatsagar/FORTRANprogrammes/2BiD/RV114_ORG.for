C RESPONSE OF 3-D ONE-STOREY BASE ISOLATED (NZ SYSTEM) STRUCTURE.
C RV114.FOR..............................................31/12/99
        DOUBLE PRECISION D1(6),V1(6),A1(6),PD(6),P1(6)
        DOUBLE PRECISION D2(6),V2(6),A2(6),PV(6),P2(6)
        DOUBLE PRECISION DD(6),DV(6),DA(6),PA(6)
        DOUBLE PRECISION EFK(6,6),EFP(6),PINC(6)
        DIMENSION X(4),Y(4),CKX(4),CKY(4),XX(6,21)
        DIMENSION SK3(3,3),SM3(3,3),CD3(3,3),OMG(3)
        DIMENSION SK6(6,6),SM6(6,6),CD6(6,6)
        DIMENSION TD(6),TV(6),TA(6),FAB(3),FD(3)
        DIMENSION WE(250),FI(3,3),FIT(3,3),PFAB(3)
        DIMENSION BS(6),PBS(6),XB(25),YB(25),SX(25),SY(25)
        DIMENSION BKX(25),BKY(25),BCX(25),BCY(25)
        DIMENSION SKB(3,3),CDB(3,3),SMB(3,3),Q(3),DQ(3,3)
        DIMENSION ZX(25),ZY(25),DZX(25),DZY(25),DZXY(25),DZYX(25)
        DIMENSION QX(25),QY(25),PQX(25),PQY(25)
        DIMENSION ALPYF(25),TK1(25),TD1(25),TK2(25),TD2(25)
        DIMENSION AA(15000)
        COMMON/RS1/SK3,SM3
        COMMON/RS2/BS,PBS,FAB,PFAB
        COMMON/RS3/EQX(6,15000),EQY(6,15000)
        COMMON/RS4/X,Y,CKX,CKY
        COMMON/RS5/SK6,SM6,CD6
        COMMON/RS6/D1,V1,A1
        COMMON/RS7/D2,V2,A2
        COMMON/RS8/PD,PV,PA
        COMMON/RS9/EFP,EFK,DD
        COMMON/RS10/TD,TV,TA
        COMMON/RS12/CD3,FM,ZETA,EX,CDABX,CDABY
        COMMON/RS13/FI,FIT
        COMMON/RS14/NB,XB,YB,BKX,BKY
        COMMON/RS15/QX,QY,PQX,PQY
        COMMON/RS16/QYF,QYD,ALP,B,T,A,NT,DT
        COMMON/RS17/AA
        OPEN(1,FILE='RV114.DAT')
        OPEN( 2,FILE='Lgpc_acc_00.txt')
        OPEN(21,FILE='Lgpc_acc_90.txt')
        OPEN(3,FILE='RV114.3')
        OPEN(4,FILE='RV114.4')
        OPEN(5,FILE='RV114.5')
C       OPEN(6,FILE='RV114.6')
C       OPEN(7,FILE='7.DAT')
C       OPEN(8,FILE='8.DAT')
C       OPEN(9,FILE='9.DAT')
C       OPEN(10,FILE='RV114.10')
        OPEN(11,FILE='RV114G.DAT')
        READ(1,*)NDT,DT,NDIV
        READ(1,*)FM,NB
        READ(1,*)N,NR
        READ(1,*)KEY,KSN,KINT,KXY
        READ(1,*)ANG,WE(1),EY
C       KEY=0,1,<SCREEN ON>,<SCREEN OFF>
C       KSN=0,1,<EARTHQUAKE>,<SINUSOIDAL>
C       KINT=0,1,<NO INTERACTION>,<INTERACTION>
C       KXY=0,1,<DIFFERENT EARTHQUAKE>,<SAME>
        DO 11 I=1,4
        READ(1,*) X(I),Y(I)
 11     CONTINUE
        DO 15 I=1,NB
        READ(1,*) XB(I),YB(I)
 15     CONTINUE
        IF(KSN.EQ.1) GO TO 104
        DO 100 KR=1,NR
 100    READ(2,*)(EQX(KR,I),I=1,NDT)
        IF(KXY.EQ.0) GO TO 102
        DO 101 KR=1,NR
 101    READ(21,*)(EQY(KR,I),I=1,NDT)
 102    DO 103 KR=1,NR
        DO 103 I=1,NDT
        EQX(KR,I)=EQX(KR,I)/100.0
 103    EQY(KR,I)=EQY(KR,I)/100.0
        NW=1
        WE(1)=0.0
        GO TO 105
 104    WRITE(3,1358)
        READ(1,*)NW,AMPX,AMPY,RFYX,PHD
        READ(1,*) (WE(I),I=1,NW)
        NR=1
 105    RSJ=0.0
        DT=DT/NDIV
        PI=ATAN(1.0)*4.0
        READ(1,*) QYD,B,T,A,NT
        DO 1050 KK=1,N
c       Because isolators are Orthotropic EXD=EYD EBXD=EBYD EFXD=EFYD
        READ(1,*)IJK,TX,ZETA,EXD,RMMB,WRWX,TB,ZETABX,EBXD,WRWXB,EFXD,QYF
        READ(11,*)GAP1,AKGAP1,CDGAP1
        GAP2=GAP1
        AKGAP2=AKGAP1
        CDGAP2=CDGAP1
        BM=FM*RMMB
        ZETABY=ZETABX
        FWWB=(2.0*PI/TB)**2.0
        WB=2.0*PI/TB
        FWWRB=FWWB*WRWXB*WRWXB
        CKAB=(FM+BM)*WB*WB
        QYF=QYF*(FM+BM)*9.81
        ALP=(QYD*CKAB)/QYF
c        ALP=1.0
        CDABX=2.0*(FM+BM)*ZETABX*WB
        CDABY=2.0*(FM+BM)*ZETABY*WB

        DO 181 I=1,NB
        CKX(I)=0.0
        CKY(I)=0.0
        BKX(I)=0.0
        BKY(I)=0.0
        ALPYF(I)=0.0
 181    CONTINUE

c       for NB=4 only
        ALPYF(1)=QYF*(1.0-ALP)*0.25*(1+2.0*EFXD)
        ALPYF(2)=QYF*(1.0-ALP)*0.25*(1-2.0*EFXD)
        ALPYF(3)=QYF*(1.0-ALP)*0.25*(1-2.0*EFXD)
        ALPYF(4)=QYF*(1.0-ALP)*0.25*(1+2.0*EFXD)

        FWWX=(2.0*PI/TX)**2.0
        AKX=FM*FWWX
        AKY=AKX
        FWWR=FWWX*WRWX*WRWX
        CALL CKXY(AKX,AKY,EX,EXD)
        CALL STIF3(SK3,AKX,AKY,EX,EY)
        CALL DAMPB(CDABX,CDABY,BCX,BCY,CDB)
        FMR=SK3(3,3)/FWWR
        ANG=ANG*PI/180.0

        CALL BKXY(CKAB,EBX,EBXD,BKX,BKY)
        SKB3=0.0
        DO 152 I=1,NB
        SKB3= SKB3 + BKX(I)*YB(I)*YB(I) + BKY(I)*XB(I)*XB(I)
  152   CONTINUE
        BMR=(SKB3/FWWRB)-FMR

        WRITE(3,1301)TX
        WRITE(3,1320)EXD
        WRITE(3,1321)RMMB
        WRITE(3,1323)WRWX
        WRITE(3,1336)FM
        WRITE(3,1337)BM
        WRITE(3,1338)FMR
        WRITE(3,1371)BMR
        WRITE(3,1303)ZETA
        WRITE(3,1304)CDABX
        WRITE(3,1305)CDABY
        WRITE(3,1372)CDB(3,3)
        WRITE(3,1339)AKX
        WRITE(3,1340)AKY
        WRITE(3,1341)SK3(3,3)
        WRITE(3,1344)EX
        WRITE(3,1345)EY
        WRITE(3,1317)
        WRITE(3,1318)
        DO 106 I=1,4
        WRITE(3,1319) I,X(I),Y(I),CKX(I),CKY(I)
 106    CONTINUE
        WRITE(3,1381)
        WRITE(3,1318)
        DO 107 I=1,NB
        WRITE(3,1319) I,XB(I),YB(I),BKX(I),BKY(I)
 107    CONTINUE
        CALL SMS3(FM,FMR,SM3)
        CALL DAMP3(FWWX,WRWX,EXD,OMG,SM3)
        CALL SMSB(BM,BMR,SMB)
        CALL SMS6(SM3,SMB,SM6)
        CALL DAMP6(CD3,CDB,CD6)
        CALL STIF6(SK3,SK6,SKB)
        WRITE(3,1325)
        WRITE(3,1302)((SM3(I,J),J=1,3),I=1,3)
        WRITE(3,1300)
        WRITE(3,1302) (OMG(I),I=1,3)
        WRITE(3,1329)
        WRITE(3,1302)((FI(I,J),J=1,3),I=1,3)
        WRITE(3,1326)
        WRITE(3,1302)((CD3(I,J),J=1,3),I=1,3)
        WRITE(3,1310)
        WRITE(3,1302)((SK3(I,J),J=1,3),I=1,3)
        WRITE(3,1385)
        WRITE(3,1302)((SMB(I,J),J=1,3),I=1,3)
        WRITE(3,1383)
        WRITE(3,1302)((CDB(I,J),J=1,3),I=1,3)
        WRITE(3,1384)
        WRITE(3,1302)((SKB(I,J),J=1,3),I=1,3)
        WRITE(3,1328)
        WRITE(3,1312)((SM6(I,J),J=1,6),I=1,6)
        WRITE(3,1327)
        WRITE(3,1312)((CD6(I,J),J=1,6),I=1,6)
        WRITE(3,1311)
        WRITE(3,1312)((SK6(I,J),J=1,6),I=1,6)
c	Vasant
        SUK4=SK6(4,4)
        SUK5=SK6(5,5)
        SUK6=SK6(6,6)
        SUD4=CD6(4,4)
        SUD5=CD6(5,5)
        SUD6=CD6(6,6)

        DO 1040 KR=1,NR
C       WRITE(*,1357) KK,KR
        DO 1030 KW=1,NW
C       WRITE(*,1359) KK,KW
        W=WE(KW)*2.0*PI/TX
C       INITIALIZATION OF VECTORS.
C       Q=CONTRIBUTION TO LOAD VECTOR
C       DQ=CONTRIBUTION TO EFFK
        DO 108 I=1,3
        Q(I)=0.0
        DO 108 J=1,3
        DQ(I,J)=0.0
 108    CONTINUE
        DO 109 I=1,NB
        ZX(I)=0.0
        ZY(I)=0.0
        DZX(I)=0.0
        DZY(I)=0.0
        DZXY(I)=0.0
        DZYX(I)=0.0
        SX(I)=1.0
        SY(I)=1.0
        QX(I)=0.0
        QY(I)=0.0
        PQY(I)=0.0
        PQX(I)=0.0
 109    CONTINUE
        DO 110 I=1,6
        D1(I)=0.0
        V1(I)=0.0
        A1(I)=0.0
        D2(I)=0.0
        V2(I)=0.0
        A2(I)=0.0
        DD(I)=0.0
        DV(I)=0.0
        DA(I)=0.0
        PD(I)=0.0
        PV(I)=0.0
        PA(I)=0.0
        P1(I)=0.0
        P2(I)=0.0
        BS(I)=0.0
        PBS(I)=0.0
        PINC(I)=0.0
 110    CONTINUE
        TIME = 0.0
        PD27=0.0
        DO 1020 K=1,NDT-1
        IF(KSN.EQ.0) GO TO 111
        EQX(KR,K)=AMPX*SIN(W*TIME)
        EQY(KR,K)=AMPY*SIN(RFYX*W*TIME + PHD)
        TT=TIME+NDIV*DT
        EQX(KR,K+1)=AMPX*SIN(W*TT)
        EQY(KR,K+1)=AMPY*SIN(RFYX*W*TT + PHD)
 111    RSJ=0.0
        P1(1)=-FM*(EQX(KR,K)*COS(ANG) - EQY(KR,K)*SIN(ANG))
        P1(2)=-FM*(EQX(KR,K)*SIN(ANG) + EQY(KR,K)*COS(ANG))
        P1(4)=-BM*(EQX(KR,K)*COS(ANG) - EQY(KR,K)*SIN(ANG))
        P1(5)=-BM*(EQX(KR,K)*SIN(ANG) + EQY(KR,K)*COS(ANG))
        P2(1)=-FM*(EQX(KR,K+1)*COS(ANG) - EQY(KR,K+1)*SIN(ANG))
        P2(2)=-FM*(EQX(KR,K+1)*SIN(ANG) + EQY(KR,K+1)*COS(ANG))
        P2(4)=-BM*(EQX(KR,K+1)*COS(ANG) - EQY(KR,K+1)*SIN(ANG))
        P2(5)=-BM*(EQX(KR,K+1)*SIN(ANG) + EQY(KR,K+1)*COS(ANG))
        DO 112 I=1,6
        PINC(I)=(P2(I)-P1(I))/NDIV
 112    CONTINUE
        DO 1010 KS=1,NDIV
        TIME = TIME + DT
        DO 113 I=1,3
        Q(I)=0.0
        DO 113 J=1,3
        DQ(I,J)=0.0
 113    CONTINUE
        DO 116 I=1,NB
        IF(KINT.EQ.1) GO TO 114
c       ALPYF=QYF*(1.0-ALP)
        VVX=V1(4)-YB(I)*V1(6)
        VVY=V1(5)+XB(I)*V1(6)
        AAX=A1(4)-YB(I)*A1(6)
        AAY=A1(5)+XB(I)*A1(6)
        AV=ABS(VVX)
        CALL WEN(SX(I),ZX(I),DZX(I),AV)
        AV=ABS(VVY)
        CALL WEN(SY(I),ZY(I),DZY(I),AV)
        X1=BKX(I)
        X2=ALPYF(I)*0.75/DT
        X3=ALPYF(I)*0.25*VVX + ALPYF(I)*0.125*AAX*DT
        X4=ALPYF(I)*0.25*VVY + ALPYF(I)*0.125*AAY*DT
        DQ(1,1)=DQ(1,1) + X1+X2*DZX(I)
        DQ(1,2)=DQ(1,2)
        DQ(1,3)=DQ(1,3) - YB(I)*(X1+X2*DZX(I))
        DQ(2,1)=DQ(2,1)
        DQ(2,2)=DQ(2,2) + X1+X2*DZY(I)
        DQ(2,3)=DQ(2,3) + XB(I)*(X1+X2*DZY(I))
        DQ(3,1)=DQ(3,1) - YB(I)*(X1+X2*DZX(I))
        DQ(3,2)=DQ(3,2) + XB(I)*(X1+X2*DZY(I))
        DQ(3,3)=DQ(3,3) + X1*(XB(I)*XB(I)+YB(I)*YB(I))
        DQ(3,3)=DQ(3,3) + X2*(DZY(I)*XB(I)*XB(I))
        DQ(3,3)=DQ(3,3) + X2*(DZX(I)*YB(I)*YB(I))
        Q(1)=Q(1) + DZX(I)*X3
        Q(2)=Q(2) + DZY(I)*X4
        Q(3)=Q(3) + XB(I)*DZY(I)*X4 - YB(I)*DZX(I)*X3
        GO TO 116
 114    RSJ=0.0
c       ALPYF=QYF*(1.0-ALP)
        VVX=V1(4)-YB(I)*V1(6)
        VVY=V1(5)+XB(I)*V1(6)
        AAX=A1(4)-YB(I)*A1(6)
        AAY=A1(5)+XB(I)*A1(6)
        AVX=ABS(VVX)
        AVY=ABS(VVY)
        CALL WENINT(SX(I),SY(I),ZX(I),DZX(I),DZXY(I),
     .  ZY(I),DZY(I),DZYX(I),AVX,AVY)
        X1=BKX(I)
        X2=ALPYF(I)*0.75/DT
        X3=ALPYF(I)*0.25*VVX + ALPYF(I)*0.125*AAX*DT
        X4=ALPYF(I)*0.25*VVY + ALPYF(I)*0.125*AAY*DT
        DQ(1,1)=DQ(1,1) + X1 + X2*DZX(I)
        DQ(1,2)=DQ(1,2) + X2*DZXY(I)
        DQ(1,3)=DQ(1,3) - YB(I)*(X1+X2*DZX(I)) + XB(I)*X2*DZXY(I)
        DQ(2,1)=DQ(2,1) + X2*DZYX(I)
        DQ(2,2)=DQ(2,2) + X1+X2*DZY(I)
        DQ(2,3)=DQ(2,3) + XB(I)*(X1+X2*DZY(I)) - YB(I)*X2*DZYX(I)
        DQ(3,1)=DQ(3,1) - YB(I)*(X1+X2*DZX(I)) + XB(I)*X2*DZXY(I)
        DQ(3,2)=DQ(3,2) + XB(I)*(X1+X2*DZY(I)) - YB(I)*X2*DZYX(I)
        DQ(3,3)=DQ(3,3) + X1*(XB(I)*XB(I)+YB(I)*YB(I))
        DQ(3,3)=DQ(3,3) + X2*(DZY(I)*XB(I)*XB(I))
        DQ(3,3)=DQ(3,3) + X2*(DZX(I)*YB(I)*YB(I))
        DQ(3,3)=DQ(3,3) - XB(I)*YB(I)*(DZXY(I)+DZYX(I))
        Q(1)=Q(1) + DZX(I)*X3 + DZXY(I)*X4
        Q(2)=Q(2) + DZY(I)*X4 + DZYX(I)*X3
        Q(3)=Q(3) + XB(I)*DZY(I)*X4 - YB(I)*DZX(I)*X3
        Q(3)=Q(3) - YB(I)*DZXY(I)*X4 + XB(I)*DZYX(I)*X3
 116    CONTINUE
        DO 117 I=1,6
        P2(I)=P1(I)+PINC(I)
 117    CONTINUE
        GAX=-P2(1)/FM
        GAY=-P2(2)/FM
        IF(K+KS-2)1020,118,120
 118    CALL ACCN(P1,FAB,FD)
        DO 119 I=1,6
        A1(I)=A2(I)
        A2(I)=0.0
 119    CONTINUE
 120    RSJ=0.0
 130    CALL EFFP(PINC,EFP,DT,Q)
        CALL EFFK(EFK,DT,DQ)
        CALL GAUS(6,EFK,EFP,DD)
 140    CALL DVBS(1,DT,DV,DD)
        DO 141 I=1,3
        FAB(I)=0.0
        FD(I)=0.0
 141    CONTINUE
        DO 142 I=1,NB
        X1=BKX(I)
c       ALPYF=QYF*(1.0-ALP)
c       FOR LRB TAKE ALPYF=0.0
c       ALPYF=0.0
        X3=ABS(V2(4)-YB(I)*V2(6))
        X4=ABS(V2(5)+XB(I)*V2(6))
        ZX(I) = ZX(I) + DZX(I)*X3 + DZXY(I)*X4
        ZY(I) = ZY(I) + DZY(I)*X4 + DZYX(I)*X3
        X3 = X1*(D2(4)-YB(I)*D2(6)) + ALPYF(I)*ZX(I)
        X4 = X1*(D2(5)+XB(I)*D2(6)) + ALPYF(I)*ZY(I)
        FAB(1) = FAB(1) + X3
        FAB(2) = FAB(2) + X4
        FAB(3) = FAB(3) - X3*YB(I) + X4*XB(I)
        QX(I) = X3
        QY(I) = X4
 142    CONTINUE
        DO 146 I=1,NB
        VX=V2(4) - YB(I)*V2(6)
        VY=V2(5) + XB(I)*V2(6)
        IF(SX(I)*VX)143,144,144
 143    SX(I)=-SX(I)
 144    IF(SY(I)*VY)145,146,146
 145    SY(I)=-SY(I)
 146    CONTINUE
c	Vasant
* Impact
        DO 195 II=1,3
        FD(II)=0.0
195     CONTINUE
        DO 196 II=1,NB
        TK1(II)=0.0
        TK2(II)=0.0
        TD1(II)=0.0
        TD2(II)=0.0
196     CONTINUE

        DO 301 I=1,NB

        XGAP=ABS(D2(4))+(10.0/2.0)*ABS(D2(6))
        YGAP=ABS(D2(5))+(10.0/2.0)*ABS(D2(6))

        IF(XGAP.GT.GAP1) THEN
        FD(1)=AKGAP1*GAP1*D2(4)/ABS(D2(4))
        FD(3)=FD(3)+(AKGAP1*GAP1*D2(4)/ABS(D2(4)))*YB(I)
        ENDIF
        IF(YGAP.GT.GAP2) THEN
        FD(2)=AKGAP2*GAP2*D2(5)/ABS(D2(5))
        FD(3)=FD(3)+(AKGAP2*GAP2*D2(5)/ABS(D2(5)))*XB(I)
        ENDIF

        CALL ACCN(P2,FAB,FD)
        CALL PEAK(TIME,P1,P2,NB,GAX,GAY)

        IF(XGAP-GAP1) 201,200,200
200     SK6(4,4)=SUK4+AKGAP1
        CD6(4,4)=SUD4+CDGAP1
        TK1(I)=AKGAP1*YB(I)*YB(I)
        TD1(I)=CDGAP1*YB(I)*YB(I)
        GO TO 202
201     SK6(4,4)=SUK4
        CD6(4,4)=SUD4
        TK1(I)=0.0
        TD1(I)=0.0

202     IF(YGAP-GAP2) 204,203,203
203     SK6(5,5)=SUK5+AKGAP2
        CD6(5,5)=SUD5+CDGAP2
        TK2(I)=AKGAP2*XB(I)*XB(I)
        TD2(I)=CDGAP2*XB(I)*XB(I)
        GO TO 205
204     SK6(5,5)=SUK5
        CD6(5,5)=SUD5
        TK2(I)=0.0
        TD2(I)=0.0

205     SK6(6,6)=SUK6+TK1(I)+TK2(I)
        CD6(6,6)=SUD6+TD1(I)+TD2(I)

301     CONTINUE

* Impact
c	Vasant
1010    CONTINUE
        IF(KEY.EQ.1) GO TO 147
        WRITE(*,1315)
        WRITE(*,1316) K,TIME,(FAB(I),I=1,3)
147     RSJ=0.0
C       WRITE(5,1354) W,(PD(I),I=1,6)
C       WRITE(5,1346) TIME,(D2(I),I=1,6)
C       WRITE(6,1346) TIME,(BS(I), I=1,3),(FAB(I),I=1,3)
C       WRITE(7,1313) D2(4),FAB(1)
C       WRITE(5,1356) IJK,(PD(I),I=1,6)

C       CORNER BEARING DISPLACEMENT
        D27=ABS(D2(5))+(10.0/2.0)*ABS(D2(6))
        IF(ABS(PD27).LT.ABS(D27)) PD27=D27
        WRITE(4,1354) TIME,AA(2)/10.0,D2(5)*100.0

c        WRITE(4,1354) TIME,A2(1)/10.0,A2(2)/10.0
c	1  ,D2(4)*100.0,D2(5)*100.0
c       WRITE(4,1354) D2(1)*100.0
c       1,D2(2)*100.0,D2(3)*100.0,D2(4)*100.0,D2(5)*100.0,D2(6)*100.0
c
c       WRITE(6,1354) A2(1)/10.0
c       1,A2(2)/10.0,A2(3)/10.0,A2(4)/10.0,A2(5)/10.0,A2(6)/10.0

c       WRITE(7,1313) D2(4)*100.0,FAB(1)
c       WRITE(8,1313) D2(5)*100.0,FAB(2)
c       WRITE(9,1313) D2(6)*100.0,FAB(3)

c       WRITE(10,1354) BS(1),BS(2),BS(3),BS(4),BS(5),BS(6)

 1020   CONTINUE
        WRITE(3,1355)
        WRITE(3,1353)IJK,TX,EXD,RMMB,WRWX,ANG
        WRITE(3,1377)
        WRITE(3,1378)QYF,QYD,ALP,B,T,A,NT
        WRITE(3,1365) PFAB(1)
        WRITE(3,1366) PFAB(2)
        WRITE(3,1373) PFAB(3)
        WRITE(3,1351) PBS(1)
        WRITE(3,1352) PBS(2)
        WRITE(3,1363) PBS(3)
        WRITE(3,1374) PBS(4)
        WRITE(3,1375) PBS(5)
        WRITE(3,1376) PBS(6)
        WRITE(3,1350)
        DO 148 I=1,6
        WRITE(3,1349) PD(I),TD(I),PV(I),TV(I),PA(I),TA(I)
 148    CONTINUE
        WRITE(3,1386)
        WRITE(3,1318)
        DO 149 I=1,NB
        WRITE(3,1319) I,XB(I),YB(I),PQX(I),PQY(I)
 149    CONTINUE
        IF(KSN.EQ.0) GO TO 150
 150    DO 151 I=1,6
        XX(I,KR)=PD(I)
 151    CONTINUE
 1030   CONTINUE
 1040   CONTINUE
        DO 156 I=1,6
        PD(I)=0.0
        DO 155 J=1,NR
        PD(I)=PD(I)+ABS(XX(I,J))
 155    CONTINUE
        PD(I)=PD(I)/NR
 156    CONTINUE
c        WRITE(5,1387) ABS(PA(1))/10.0,ABS(PD(4))*100.0
c        WRITE(5,1387) PD(4)*100.0,PD(5)*100.0,PD(6)*100.0,PD27*100.0
c        WRITE(5,1387) ABS(PA(1))/10.0,ABS(PA(2))/10.0,ABS(PA(3))/10.0
c        WRITE(5,1388) ABS(PD(4))*100.0,ABS(PD(5))*100.0,ABS(PD(6))*100.0

        WRITE(5,1387) ABS(PA(1))/10.0,ABS(PA(2))/10.0,ABS(PA(3))/10.0
        WRITE(5,1387) ABS(PD(4))*100.0,ABS(PD(5))*100.0,ABS(PD(6))*100.0
        WRITE(5,1387) PD27*100.0
c       WRITE(5,1356) IJK,(PA(I)/10.0,I=1,6)
c       WRITE(5,1356) IJK,(PD(I)*100.0,I=1,6)
 1050   CONTINUE
C
 1300   FORMAT(/,10X,'FIXED BASE NATURAL FREQUENCIES ')
 1301   FORMAT(3X,'FIXED BASE UNCOUPLED TIME PERIOD = ',E13.5)
 1302   FORMAT(5X,3E13.5)
 1303   FORMAT(3X,'SUPERSTRUCTURE DAMPING CONSTANT  = ',E13.5)
 1304   FORMAT(3X,'DAMPING OF ISOLATOR IN       X   = ',E13.5)
 1305   FORMAT(3X,'DAMPING OF ISOLATOR IN       Y   = ',E13.5)
 1306   FORMAT(/,10X,'FIXED BASE MODAL COLUMN MATRIX')
 1307   FORMAT(/,10X,'FIXED BASE STIFFNESS MATRIX')
 1309   FORMAT(/,10X,'INVERSE ASSEMBLED MASS MATRIX')
 1310   FORMAT(/,10X,'FIXED BASE STIFFNESS MATRIX')
 1311   FORMAT(/,10X,'ASSEMBLED STIFFNESS MATRIX')
 1312   FORMAT(5X,6E13.5)
 1313   FORMAT(3X,4E13.5)
 1315   FORMAT(6X,'NDT   TIME       FABX  ',8X,'FABY',11X,'FABZ')
 1316   FORMAT(5X,I4,F8.3,5X,3E13.5)
 1317   FORMAT(/,3X,'C.N.',3X,'X-CORD',3X,'Y-CORD',6X,'CKX',6X,'CKY')
 1318   FORMAT(3X,'----',3X,'------',3X,'------',6X,'---',6X,'---')
 1319   FORMAT(3X,I2,2F9.2,2X,2F10.2)
 1320   FORMAT(3X,'ECEENTRICITY RATIO  X            = ',E13.5)
 1321   FORMAT(3X,'MASS RATIO FLOOR TO BASE         = ',E13.5)
 1322   FORMAT(3X,'BASE TO LATERAL FREQUENCY RATIO  = ',E13.5)
 1323   FORMAT(3X,'TORSIONAL TO LATERAL FREQ. RATIO = ',E13.5)
 1324   FORMAT(3X,'ECEENTRICITY RATIO  Y            = ',E13.5)
 1325   FORMAT(/,10X,'FIXED BASE MASS MATRIX')
 1326   FORMAT(/,10X,'FIXED BASE DAMPING MATRIX')
 1327   FORMAT(/,10X,'ASSEMBLED DAMPING MATRIX')
 1328   FORMAT(/,10X,'ASSEMBLED MASS MATRIX')
 1329   FORMAT(/,10X,'FIXED BASE MODE-SHAPE MATRIX')
 1336   FORMAT(3X,'MASS OF FLOOR                    = ',E13.5)
 1337   FORMAT(3X,'MASS OF BASE                     = ',E13.5)
 1338   FORMAT(3X,'ROTATIONAL MASS OF FLOOR         = ',E13.5)
 1339   FORMAT(3X,'STIFFNESS OF COLUMN IN X         = ',E13.5)
 1340   FORMAT(3X,'STIFFNESS OF COLUMN IN Y         = ',E13.5)
 1341   FORMAT(3X,'ROTATIONAL STIFFNESS OF FLOOR    = ',E13.5)
 1342   FORMAT(3X,'TOTAL STIFFNESS OF ISOLATOR IN X = ',E13.5)
 1343   FORMAT(3X,'TOTAL STIFFNESS OF ISOLATOR IN Y = ',E13.5)
 1344   FORMAT(3X,'ECCENTRICITY IN X DIRECTION      = ',E13.5)
 1345   FORMAT(3X,'ECCENTRICITY IN Y DIRECTION      = ',E13.5)
 1346   FORMAT(F6.3,1X,7E10.3)
 1349   FORMAT(3X,E10.4,2X,F6.2,E10.4,2X,F6.2,E10.4,2X,F6.2)
 1350   FORMAT(7X,'PEAK DISP  T     PEAK VEL.   T    PEAK ACCN    T')
 1351   FORMAT(10X,'MAXIMUM BASE SHEAR     X     =  ',E13.5)
 1352   FORMAT(10X,'MAXIMUM BASE SHEAR     Y     =  ',E13.5)
 1353   FORMAT(I3,6X,6F6.2,F6.2,F6.2)
 1354   FORMAT(6E15.4)
 1355   FORMAT(11X,'TX   EX/D  M/MB  WR/WX  ANG')
 1356   FORMAT(I10,6E14.6)
 1357   FORMAT(10X,'N = ',I3,4X,'NR = ',I3)
 1358   FORMAT(10X,'RESPONSE TO SINOSOIDAL EXCITATION ',/)
 1359   FORMAT(10X,'N = ',I3,4X,'NW = ',I3)
 1361   FORMAT(/,10X,'A  MATRIX ')
 1362   FORMAT(3X,6F7.2,4F8.2,2F5.2)
 1363   FORMAT(10X,'MAXIMUM BASE TORQUE    Z     =  ',E13.5)
 1365   FORMAT(10X,'MAXIMUM FORCE IN  ISOLATOR X =  ',E13.5)
 1366   FORMAT(10X,'MAXIMUM FORCE IN  ISOLATOR Y =  ',E13.5)
 1367   FORMAT(I4,F10.3,6E14.6)
 1371   FORMAT(3X,'BASE MASS ROTATIONAL             = ',E13.5)
 1372   FORMAT(3X,'DAMPING OF ISOLATOR IN       Z   = ',E13.5)
 1373   FORMAT(10X,'MAXIMUM FORCE IN  ISOLATOR Z =  ',E13.5)
 1374   FORMAT(10X,'MAXIMUM BASE SHEAR     BX    =  ',E13.5)
 1375   FORMAT(10X,'MAXIMUM BASE SHEAR     BY    =  ',E13.5)
 1376   FORMAT(10X,'MAXIMUM BASE TORQUE    BZ    =  ',E13.5)
 1377   FORMAT(11X,'QYF  QYD   ALP    B       T      A     N ')
 1378   FORMAT(9X,6F6.2,I6)
 1380   FORMAT(3X,'DECK TO BASE TORSNL. FREQ. RATIO = ',E13.5)
 1381   FORMAT(/,3X,'B.N.',3X,'X-CORD',3X,'Y-CORD',6X,'BKX',6X,'BKY')
 1382   FORMAT(3X,'TOTAL STIFFNESS OF ISOLATOR IN Z = ',E13.5)
 1383   FORMAT(/,10X,'BASE DAMPING MATRIX')
 1384   FORMAT(/,10X,'BASE STIFFNESS MATRIX')
 1385   FORMAT(/,10X,'BASE MASS MATRIX')
 1386   FORMAT(/,3X,'B.N.',3X,'X-CORD',3X,'Y-CORD',6X,'PQX',6X,'PQY')
 1387   FORMAT(6F10.4)
c 1387   FORMAT(F8.3,F8.3,F9.4)
c 1388   FORMAT(F8.2,F8.2,F9.3)
        STOP
        END
 
        SUBROUTINE CKXY(AKX,AKY,EX,EXD)
        DIMENSION X(4),Y(4),CKX(4),CKY(4)
        COMMON/RS4/X,Y,CKX,CKY
        D=X(1)-X(2)
        EX=EXD*D
        CKX(1)=0.25*AKX*(1.0 + 2.0*EXD)
        CKX(2)=0.25*AKX*(1.0 - 2.0*EXD)
        CKX(3)=0.25*AKX*(1.0 - 2.0*EXD)
        CKX(4)=0.25*AKX*(1.0 + 2.0*EXD)
        DO 1 I=1,4
        CKY(I)=CKX(I)
   1    CONTINUE
        RETURN
        END
 
        SUBROUTINE STIF3(SK3,AKX,AKY,EX,EY)
        DIMENSION X(4),Y(4),CKX(4),CKY(4)
        DIMENSION SK3(3,3)
        COMMON/RS4/X,Y,CKX,CKY
        DO 1 I=1,3
        DO 1 J=1,3
        SK3(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,4
        SK3(1,1) = SK3(1,1) + CKX(I)
        SK3(1,2) = SK3(1,2)
        SK3(1,3) = SK3(1,3) - Y(I)*CKX(I)
        SK3(2,1) = SK3(2,1)
        SK3(2,2) = SK3(2,2) + CKY(I)
        SK3(2,3) = SK3(2,3) + X(I)*CKY(I)
        SK3(3,1) = SK3(3,1) - Y(I)*CKX(I)
        SK3(3,2) = SK3(3,2) + X(I)*CKY(I)
        SK3(3,3) = SK3(3,3) + CKX(I)*Y(I)*Y(I) + CKY(I)*X(I)*X(I)
   2    CONTINUE
C       WRITE(3,1310)
C       WRITE(3,1302)((SK3(I,J),J=1,3),I=1,3)
C 1310  FORMAT(/,10X,'FIXED BASE STIFFNESS MATRIX')
C 1302  FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE SMS3(FM,FMR,SM3)
        DIMENSION SM3(3,3)
        DO 1 I=1,3
        DO 1 J=1,3
        SM3(I,J)=0.0
   1    CONTINUE
        SM3(1,1) = SM3(1,1) + FM
        SM3(2,2) = SM3(2,2) + FM
        SM3(3,3) = SM3(3,3) + FMR
C       WRITE(3,1325)
C       WRITE(3,1302)((SM3(I,J),J=1,3),I=1,3)
C 1325  FORMAT(/,10X,'FIXED BASE MASS MATRIX')
C 1302  FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE DAMP3(WW0,ALPHA,BETA,OMG,SM3)
        DIMENSION OMG(3),GM(3),CD3(3,3),FI(3,3)
        DIMENSION FIT(3,3),SM3(3,3)
        COMMON/RS12/CD3,FM,ZETA,EX,CDABX,CDABY
        COMMON/RS13/FI,FIT
        IF(EX-0.00001) 10,10,1
   1    OMG(1)=WW0
        A1 = 1.0 + ALPHA*ALPHA
        A2 = 1.0 - ALPHA*ALPHA
        A3 = ALPHA*ALPHA*BETA*BETA
        A4 = A2*A2 + 8.0*A3
        OMG(2) = WW0*(A1 - SQRT(A4))/2.0
        OMG(3) = WW0*(A1 + SQRT(A4))/2.0
        A5 = (1.0 - OMG(2)/WW0)/EX
        A6 = (1.0 - OMG(3)/WW0)/EX
        A7 = (A6-A5)**2
        A8 = (EX/(BETA*ALPHA))**2
        GM(1)=FM
        GM(2)=FM*(1.0 + A8*A5*A5/2.0)
        GM(3)=FM*(1.0 + A8*A6*A6/2.0)
        DO 2 I=1,3
        OMG(I)=SQRT(OMG(I))
        DO 2 J=1,3
        CD3(I,J) = 0.0
        FI(I,J) = 0.0
        FIT(I,J) = 0.0
   2    CONTINUE
        FI(1,1) = FI(1,1) + 1.0
        FI(2,2) = FI(2,2) + 1.0
        FI(2,3) = FI(2,3) + 1.0
        FI(3,2) = FI(3,2) + A5
        FI(3,3) = FI(3,3) + A6
        CD3(1,1) = CD3(1,1) + GM(1)*OMG(1)*A7
        A10 = GM(2)*OMG(2)*A6
        A11 = GM(3)*OMG(3)*A5
        CD3(2,2) = CD3(2,2) + A10*A6 + A11*A5
        CD3(2,3) = CD3(2,3) + A10 + A11
        CD3(3,3) = CD3(3,3) + A10/A6 + A11/A5
        CD3(3,2) = CD3(3,2) + CD3(2,3)
        DO 3 I=1,3
        DO 3 J=1,3
        CD3(I,J) = 2.0*ZETA*CD3(I,J)/A7
   3    CONTINUE
        GO TO 13
  10    DO 4 I=1,3
        OMG(I)=0.0
   4    CONTINUE
        OMG(1)=WW0
        OMG(2)=WW0
        OMG(3)=WW0*ALPHA*ALPHA
        DO 11 I=1,3
        OMG(I)=SQRT(OMG(I))
        DO 11 J=1,3
        CD3(I,J) = 0.0
        FI(I,J) = 0.0
        FIT(I,J) = 0.0
  11    CONTINUE
        FI(1,1) = FI(1,1) + 1.0
        FI(2,2) = FI(2,2) + 1.0
        FI(3,3) = FI(3,3) + 1.0
        DO 12 J=1,3
        CD3(J,J) = CD3(J,J) + 2.0*ZETA*SM3(J,J)*OMG(J)
  12    CONTINUE
  13    DO 14 I=1,3
        DO 14 J=1,3
        FIT(I,J)=FI(J,I)
  14    CONTINUE
C       WRITE(3,1300)
C       WRITE(3,1302) (OMG(I),I=1,3)
C       WRITE(3,1329)
C       WRITE(3,1302)((FI(I,J),J=1,3),I=1,3)
C       WRITE(3,1326)
C       WRITE(3,1302)((CD3(I,J),J=1,3),I=1,3)
C 1300  FORMAT(/,10X,'FIXED BASE NATURAL FREQUENCIES ')
C 1326  FORMAT(/,10X,'FIXED BASE DAMPING MATRIX')
C 1329  FORMAT(/,10X,'FIXED BASE MODE-SHAPE MATRIX')
C 1302  FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE SMSB(BM,BMR,SMB)
        DIMENSION SMB(3,3)
        DO 1 I=1,3
        DO 1 J=1,3
        SMB(I,J)=0.0
   1    CONTINUE
        SMB(1,1) = SMB(1,1) + BM
        SMB(2,2) = SMB(2,2) + BM
        SMB(3,3) = SMB(3,3) + BMR
C       WRITE(3,1385)
C       WRITE(3,1302)((SMB(I,J),J=1,3),I=1,3)
C 1385  FORMAT(/,10X,'BASE MASS MATRIX')
C 1302  FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE DAMPB(CDABX,CDABY,BCX,BCY,CDB)
        DIMENSION XB(25),YB(25),BKX(25),BKY(25)
        DIMENSION CDB(3,3),BCX(25),BCY(25)
        COMMON/RS14/NB,XB,YB,BKX,BKY
        DO 1 I=1,NB
        BCX(I)=CDABX/NB
        BCY(I)=CDABY/NB
   1    CONTINUE
        DO 2 I=1,3
        DO 2 J=1,3
        CDB(I,J)=0.0
   2    CONTINUE
        DO 3 I=1,NB
        CDB(1,1) = CDB(1,1) + BCX(I)
        CDB(1,2) = CDB(1,2)
        CDB(1,3) = CDB(1,3) - YB(I)*BCX(I)
        CDB(2,1) = CDB(2,1)
        CDB(2,2) = CDB(2,2) + BCY(I)
        CDB(2,3) = CDB(2,3) + XB(I)*BCY(I)
        CDB(3,1) = CDB(3,1) - YB(I)*BCX(I)
        CDB(3,2) = CDB(3,2) + XB(I)*BCY(I)
        CDB(3,3) = CDB(3,3) + BCX(I)*YB(I)*YB(I) + BCY(I)*XB(I)*XB(I)
   3    CONTINUE
C       WRITE(3,1383)
C       WRITE(3,1302)((CDB(I,J),J=1,3),I=1,3)
C 1383  FORMAT(/,10X,'BASE DAMPING MATRIX')
C 1302  FORMAT(5X,3E13.5)
        RETURN
        END

        SUBROUTINE BKXY(CKAB,EBX,EBXD,BKX,BKY)
        DIMENSION XB(25),BKX(25),BKY(25)
        COMMON/RS14/NB,XB,YB
        DB=XB(1)-XB(2)
        EBX=EBXD*DB
c ****        WARNING: for NB=4 only
        BKX(1)=0.25*CKAB*(1.0 + 2.0*EBXD)
        BKX(2)=0.25*CKAB*(1.0 - 2.0*EBXD)
        BKX(3)=0.25*CKAB*(1.0 - 2.0*EBXD)
        BKX(4)=0.25*CKAB*(1.0 + 2.0*EBXD)
        DO 1 I=1,4
1       BKY(I)=BKX(I)
        RETURN
        END

        SUBROUTINE STIFB(SKB)
        DIMENSION XB(25),YB(25),BKX(25),BKY(25)
        DIMENSION SKB(3,3)
        COMMON/RS14/NB,XB,YB,BKX,BKY
        DO 1 I=1,3
        DO 1 J=1,3
        SKB(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,NB
        SKB(1,1) = SKB(1,1) + BKX(I)
        SKB(1,2) = SKB(1,2)
        SKB(1,3) = SKB(1,3) - YB(I)*BKX(I)
        SKB(2,1) = SKB(2,1)
        SKB(2,2) = SKB(2,2) + BKY(I)
        SKB(2,3) = SKB(2,3) + XB(I)*BKY(I)
        SKB(3,1) = SKB(3,1) - YB(I)*BKX(I)
        SKB(3,2) = SKB(3,2) + XB(I)*BKY(I)
        SKB(3,3) = SKB(3,3) + BKX(I)*YB(I)*YB(I) + BKY(I)*XB(I)*XB(I)
   2    CONTINUE
C       WRITE(3,1384)
C       WRITE(*,1302)((SKB(I,J),J=1,3),I=1,3)
C1384   FORMAT(/,10X,'BASE STIFFNESS MATRIX')
C1302   FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE DAMP6(CD3,CDB,CD6)
        DIMENSION CD6(6,6),CD3(3,3),CDB(3,3)
        DO 1 I=1,6
        DO 1 J=1,6
        CD6(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,3
        DO 2 J=1,3
        CD6(I,J) = CD3(I,J)
   2    CONTINUE
        DO 3 I=1,3
        II=I+3
        DO 3 J=1,3
        CD6(II,J) = -CD3(I,J)
   3    CONTINUE
        DO 4 I=1,3
        II=I+3
        DO 4 J=1,3
        JJ=J+3
        CD6(II,JJ) = CDB(I,J)
   4    CONTINUE
C       WRITE(3,1327)
C       WRITE(3,1312)((CD6(I,J),J=1,6),I=1,6)
C 1327  FORMAT(/,10X,'ASSEMBLED DAMPING MATRIX')
C 1312  FORMAT(5X,6E13.5)
        RETURN
        END
 
        SUBROUTINE SMS6(SM3,SMB,SM6)
        DIMENSION SM6(6,6),SM3(3,3),SMB(3,3)
        DO 1 I=1,6
        DO 1 J=1,6
        SM6(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,3
        DO 2 J=1,3
        SM6(I,J) = SM3(I,J)
   2    CONTINUE
        DO 3 I=1,3
        DO 3 J=1,3
        JJ=J+3
        SM6(I,JJ) = SM3(I,J)
   3    CONTINUE
        DO 4 I=1,3
        II=I+3
        DO 4 J=1,3
        JJ=J+3
        SM6(II,JJ) = SMB(I,J)
   4    CONTINUE
C       WRITE(3,1328)
C       WRITE(3,1312)((SM6(I,J),J=1,6),I=1,6)
C 1328  FORMAT(/,10X,'ASSEMBLED MASS MATRIX')
C 1312  FORMAT(5X,6E13.5)
        RETURN
        END
 
        SUBROUTINE STIF6(SK3,SK6,SKB)
        DIMENSION SK6(6,6),SK3(3,3),SKB(3,3)
        DO 1 I=1,6
        DO 1 J=1,6
        SK6(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,3
        DO 2 J=1,3
        SK6(I,J) = SK3(I,J)
   2    CONTINUE
        DO 3 I=1,3
        II=I+3
        DO 3 J=1,3
        SK6(II,J) = -SK3(I,J)
   3    CONTINUE
C       WRITE(3,1311)
C       WRITE(*,1312)((SK6(I,J),J=1,6),I=1,6)
C 1311  FORMAT(/,10X,'ASSEMBLED STIFFNESS MATRIX')
C 1312  FORMAT(5X,6E13.5)
        RETURN
        END
 
        SUBROUTINE EFFP(PINC,EFP,T,Q)
        DOUBLE PRECISION EFP(6),PINC(6),D1(6),V1(6),A1(6)
        DIMENSION SK6(6,6),SM6(6,6),CD6(6,6),Q(3)
        COMMON/RS5/SK6,SM6,CD6
        COMMON/RS6/D1,V1,A1
        DO 1 I=1,6
        EFP(I) =0.0
   1    CONTINUE
        DO 2 I=1,6
        EFP(I)=EFP(I)+PINC(I)
        DO 3 J=1,6
        X1=SM6(I,J)*(6.0*V1(J)/T + 3.0*A1(J))
        X2=CD6(I,J)*(3.0*V1(J) + 0.5*T*A1(J))
        EFP(I)=EFP(I) + X1 + X2
   3    CONTINUE
   2    CONTINUE
        DO 4 I=1,3
        II=I+3
        EFP(II)=EFP(II)-Q(I)
   4    CONTINUE
C       WRITE(8,11)
C       WRITE(8,12)(EFP(I),I=1,3)
C 11    FORMAT(10X,'EFFECTIVE FORCE')
        RETURN
        END

        SUBROUTINE EFFK(EFK,T,DQ)
        DOUBLE PRECISION EFK(6,6)
        DIMENSION SK6(6,6),SM6(6,6),CD6(6,6),DQ(3,3)
        COMMON/RS5/SK6,SM6,CD6
        DO 1 I=1,6
        DO 1 J=1,6
   1    EFK(I,J)=SK6(I,J)
        DO 2 I=1,6
        DO 2 J=1,6
        EFK(I,J)=EFK(I,J) + CD6(I,J)*3.0/T + SM6(I,J)*6.0/(T*T)
   2    CONTINUE
        DO 3 I=1,3
        II=I+3
        DO 3 J=1,3
        JJ=J+3
        EFK(II,JJ)=EFK(II,JJ)+DQ(I,J)
   3    CONTINUE
C       WRITE(8,11)
C       WRITE(8,12)((EFK(I,J),J=1,3),I=1,3)
C 11    FORMAT(10X,'EFFECTIVE STIFFNESS MATRIX')
C 12    FORMAT(5X,3E13.5)
        RETURN
        END
 
        SUBROUTINE GAUS(N,A,F,X)
        DOUBLE PRECISION A(6,6),F(6),X(6)
C       CHOOSE THE K TH UNKNOWN TO BE
C       ELIMINATED FROM ALL THE SUCCEEDING EQUATION
        DO 10 K=1,(N-1)
C       ELIMINATE THE KTH UNKNOWN FROM ALL
C       (K+1) TO N EQUATIONS
        DO 20 I=(K+1),N
        QUOT=A(I,K)/A(K,K)
C       MODIFY THE LOAD VECTOR BY
        F(I)=F(I)-QUOT*F(K)
C       MODIFY THE I TH EQUATION ELEMENTS
        DO 30 J=(K+1),N
        A(I,J)=A(I,J)-QUOT*A(K,J)
  30    CONTINUE
  20    CONTINUE
  10    CONTINUE
        X(N)=F(N)/A(N,N)
        DO 40 I=(N-1),1,-1
        SUM=0.0
        DO 50 J=(I+1),N
        SUM=SUM+A(I,J)*X(J)
  50    CONTINUE
        X(I)=(F(I)-SUM)/A(I,I)
  40    CONTINUE
C       WRITE(8,11)
C       WRITE(8,12)(X(I),I=1,3)
C 11    FORMAT(10X,'DELTA DISPLACEMENT')
C 12    FORMAT(3X,3E13.5)
        RETURN
        END
 
        SUBROUTINE ACCN(P,FAB,FD)
        DOUBLE PRECISION D2(6),V2(6),A2(6),P(6)
        DIMENSION SK6(6,6),SM6(6,6),CD6(6,6)
        DIMENSION FAB(3),FC(6),FK(6),EP1(6),FD(3)
        COMMON/RS5/SK6,SM6,CD6
        COMMON/RS7/D2,V2,A2
        DO 2 I=1,6
        FC(I)=0.0
        FK(I)=0.0
        DO 1 J=1,6
        FC(I)=FC(I)+CD6(I,J)*V2(J)
        FK(I)=FK(I)+SK6(I,J)*D2(J)
   1    CONTINUE
        EP1(I)=P(I)-FC(I)-FK(I)
   2    CONTINUE
        DO 3 I=1,3
        II=I+3
c	Vasant
        EP1(II)=EP1(II)-FAB(I)+FD(I)
   3    CONTINUE
        DO 4 I=4,6
        A2(I)=EP1(I)/SM6(I,I)
   4    CONTINUE
        A2(1)=(EP1(1)-SM6(1,4)*A2(4))/SM6(1,1)
        A2(2)=(EP1(2)-SM6(2,5)*A2(5))/SM6(2,2)
        A2(3)=(EP1(3)-SM6(3,6)*A2(6))/SM6(3,3)
        RETURN
        END
 
        SUBROUTINE DVBS(KSP,DT,DV,DD)
        DOUBLE PRECISION D1(6),V1(6),A1(6)
        DOUBLE PRECISION D2(6),V2(6),A2(6)
        DOUBLE PRECISION DD(6),DV(6)
        DIMENSION SK6(6,6),SM6(6,6),CD6(6,6)
        DIMENSION BS(6),PBS(6),FAB(3),PFAB(3)
        COMMON/RS2/BS,PBS,FAB,PFAB
        COMMON/RS5/SK6,SM6,CD6
        COMMON/RS6/D1,V1,A1
        COMMON/RS7/D2,V2,A2
        DO 2 I=1,6
        D2(I)=D1(I)+DD(I)
        IF(KSP.EQ.0) GO TO 1
        DV(I)=3.0*DD(I)/DT-3.0*V1(I)-DT*A1(I)*0.5
   1    V2(I)=V1(I)+DV(I)
   2    CONTINUE
        DO 3 I=1,3
        II=I+3
        BS(I) = 0.0
        BS(II) = FAB(I)
        DO 3 J=1,3
        JJ=J+3
        BS(I) = BS(I) + SK6(I,J)*D2(J) + CD6(I,J)*V2(J)
        BS(II) = BS(II) + CD6(II,JJ)*V1(JJ)
   3    CONTINUE
        RETURN
        END
 
C       TO INVERSE A GIVEN MATRIX.
        SUBROUTINE INV(A,B,N)
        DIMENSION A(N,N),B(N,N)
        DO 1 I=1,N
   1    B(I,I)=1.0
        DO 6 I=1,N
        TEMP=A(I,I)
        DO 2 J=1,N
        A(I,J)=A(I,J)/TEMP
        B(I,J)=B(I,J)/TEMP
   2    CONTINUE
        DO 5 J=1,N
        IF(I-J) 3,5,3
   3    RT=A(J,I)
        DO 4 K=1,N
        A(J,K)=A(J,K)-A(I,K)*RT
   4    B(J,K)=B(J,K)-B(I,K)*RT
   5    CONTINUE
   6    CONTINUE
        RETURN
        END
 
        SUBROUTINE PEAK(T,P1,P2,NB,EQX,EQY)
        DOUBLE PRECISION D1(6),V1(6),A1(6),PD(6)
        DOUBLE PRECISION D2(6),V2(6),A2(6),PV(6)
        DOUBLE PRECISION P1(6),P2(6),PA(6)
        DIMENSION TD(6),TV(6),TA(6),AA(6)
        DIMENSION BS(6),PBS(6),FAB(3),PFAB(3)
        DIMENSION QX(25),QY(25),PQX(25),PQY(25)
        COMMON/RS6/D1,V1,A1
        COMMON/RS7/D2,V2,A2
        COMMON/RS8/PD,PV,PA
        COMMON/RS10/TD,TV,TA
        COMMON/RS2/BS,PBS,FAB,PFAB
        COMMON/RS15/QX,QY,PQX,PQY
        COMMON/RS17/AA
        AA(1)=A2(1)+A2(4)+EQX
        AA(2)=A2(2)+A2(5)+EQY
        AA(3)=A2(3)+A2(6)
        AA(4)=A2(4)+EQX
        AA(5)=A2(5)+EQY
        AA(6)=A2(6)
        DO 4 I=1,3
        IF(ABS(PFAB(I))-ABS(FAB(I)))3,4,4
   3    PFAB(I)=FAB(I)
   4    CONTINUE
   8    DO 10 I=1,6
        IF(ABS(PD(I))-ABS(D2(I)))9,10,10
   9    PD(I)=D2(I)
        TD(I)=T
  10    CONTINUE
        DO 12 I=1,6
        IF(ABS(PV(I))-ABS(V2(I)))11,12,12
  11    PV(I)=V2(I)
        TV(I)=T
  12    CONTINUE
        DO 14 I=1,6
        IF(ABS(PA(I))-ABS(AA(I)))13,14,14
  13    PA(I)=AA(I)
        TA(I)=T
  14    CONTINUE
        DO 16 I=1,6
        IF(ABS(PBS(I))-ABS(BS(I)))15,16,16
  15    PBS(I)=BS(I)
  16    CONTINUE
        DO 17 I=1,6
        A1(I)=A2(I)
        V1(I)=V2(I)
        D1(I)=D2(I)
        P1(I)=P2(I)
  17    CONTINUE
        DO 24 I=1,NB
        IF(ABS(PQX(I))-ABS(QX(I)))21,22,22
  21    PQX(I)=QX(I)
  22    IF(ABS(PQY(I))-ABS(QY(I)))23,24,24
  23    PQY(I)=QY(I)
  24    CONTINUE
        RETURN
        END

        SUBROUTINE WEN(V,Z,DZ,AV)
        COMMON/RS16/FY,G1,ALP,BT,G2,A,NT,DT
        FDT=DT/G1
        FZ=Z
        B1=FDT*(-G2*ABS(V)*FZ*(ABS(FZ)**(NT-1)))
     .  +  FDT*(-BT*V*(ABS(FZ)**NT) + A*V)
        FZ=Z+B1*AV/3.0
        B2=FDT*(-G2*ABS(V)*FZ*(ABS(FZ)**(NT-1)))
     .  +  FDT*(-BT*V*(ABS(FZ)**NT) + A*V)
        FZ=Z+2.0*AV*B2/3
        B3=FDT*(-G2*ABS(V)*FZ*(ABS(FZ)**(NT-1)))
     .  +  FDT*(-BT*V*(ABS(FZ)**NT) + A*V)
        DZ=(B1+3.0*B3)/4.0
        RETURN
        END

        SUBROUTINE WENINT(VX,VY,ZX,DZX,DZXY,ZY,DZY,DZYX,AVX,AVY)
        COMMON/RS16/FY,G1,ALP,BT,G2,A,NT,DT
        FDT=DT/G1
        FZX=ZX
        FZY=ZY
        A11=FDT*(-G2*ABS(VX*FZX)*FZX - BT*VX*FZX*FZX + A*VX)
        A12=FDT*(-G2*ABS(VY*FZY)*FZX - BT*VY*FZX*FZY)
        B11=FDT*(-G2*ABS(VY*FZY)*FZY - BT*VY*FZY*FZY + A*VY)
        B12=FDT*(-G2*ABS(VX*FZX)*FZY - BT*VX*FZX*FZY)
        FZX=ZX + AVX*(A11+A12)/3.0
        FZY=ZY + AVY*(B11+B12)/3.0
        A21=FDT*(-G2*ABS(VX*FZX)*FZX - BT*VX*FZX*FZX + A*VX)
        A22=FDT*(-G2*ABS(VY*FZY)*FZX - BT*VY*FZX*FZY)
        B21=FDT*(-G2*ABS(VY*FZY)*FZY - BT*VY*FZY*FZY + A*VY)
        B22=FDT*(-G2*ABS(VX*FZX)*FZY - BT*VX*FZX*FZY)
        FZX=ZX + AVX*2.0*(A21+A22)/3.0
        FZY=ZY + AVY*2.0*(B21+B22)/3.0
        A31=FDT*(-G2*ABS(VX*FZX)*FZX - BT*VX*FZX*FZX + A*VX)
        A32=FDT*(-G2*ABS(VY*FZY)*FZX - BT*VY*FZX*FZY)
        B31=FDT*(-G2*ABS(VY*FZY)*FZY - BT*VY*FZY*FZY + A*VY)
        B32=FDT*(-G2*ABS(VX*FZX)*FZY - BT*VX*FZX*FZY)
        DZX=(A11+3.0*A31)/4.0
        DZY=(B11+3.0*B31)/4.0
        DZXY=(A12+3.0*A32)/4.0
        DZYX=(B12+3.0*B32)/4.0
        RETURN
        END
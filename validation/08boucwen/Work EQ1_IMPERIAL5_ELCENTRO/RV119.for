C RESPONSE OF MULTI-STOREY BUILDING WITH N-Z SYSTEMS.
C RV119.FOR..................................31/12/02
        DIMENSION D1(11),V1(11),A1(11),PD(11),PA(11)
        DIMENSION D2(11),V2(11),A2(11),AA(11),PV(11)
        DIMENSION SK(11,11),SM(11,11),CD(11,11)
        DIMENSION DD(11),DV(11),P1(11),P2(11),PINC(11)
        DIMENSION EFK(11,11),EFP(11),WW(50)
        DIMENSION V(11,11),VI(11,11),W(11)
        DIMENSION AM(11),AK(11),ZETA(11)
        COMMON/RS1/XG(13000)
        COMMON/RS5/SM,SK,CD
        COMMON/RS6/D1,V1,A1
        COMMON/RS7/D2,V2,A2
        COMMON/RS8/PD,PV,PA
        COMMON/RS12/FY,G1,ALP,BT,G2,A,NT,DT
	   
        OPEN(1,FILE='RV119.DAT')

           OPEN(2,FILE='Imperial5_Elcentro.txt')
C   	   OPEN(2,FILE='a3Sin0p3t.txt')
C   	   OPEN(2,FILE='Landers_Lucerne.txt')
C   	   OPEN(2,FILE='Northridge_Newhall.txt')
C  	   OPEN(2,FILE='Northridge_Rinaldi.txt')
C         OPEN(2,FILE='Northridge_Sylmar.txt')

        OPEN(3,FILE='RV119.3')
        OPEN(4,FILE='RV119.41')
        OPEN(5,FILE='RV119.5')
        OPEN(6,FILE='RV119.6')
        OPEN(7,FILE='TEST1.txt')
        OPEN(8,FILE='DEBUG.txt')

        READ(1,*) NDT,DTN,NDIV,NIT
        READ(1,*) NST
        READ(1,*) (AM(I),I=1,NST)
        READ(1,*) (AK(I),I=1,NST)
        READ(1,*) (ZETA(I),I=1,NST)
        READ(1,*) BT,G2,A,NT 
        READ(1,*) KNOR,KSN
        READ(1,*) N

        IF(KSN.EQ.1) READ(1,*) AMP,NW,(WW(I),I=1,NW)
        IF(KSN.EQ.0) READ(2,*) (XG(I),I=1,NDT)
        IF(KSN.EQ.0) NW=1 

	   DO 1 I=1,NDT
1	   XG(I)=XG(I)/100.0

        PI=ATAN(1.0)*4.0

        DO 1050 KN=1,N 
        READ(1,*) IJK,RMBM,BEFF,NST,T1,F0,TB,G1

	   WRITE(*,*)IJK
        BM=RMBM*AM(1)
        W1=2.0*PI/T1
        NDOF=NST+1
        TM=BM
        DO 106 I=1,NST
        TM=TM+AM(I)
 106    CONTINUE
        WB=2.0*PI/TB
	   CKAB=TM*WB*WB
	  
	 
	   FY = F0 * TM * 9.81 
C       Q=(BEFF*PI*CKAB*D*D)/(2.0*(D-G1))

C       CKAB=(CKAB*D-Q)/D
	   Q = FY- CKAB*G1
C       FY = Q + CKAB*G1
       
C	   Q=(FY-CKAB*G1)*(FY-CKAB*G1)/(FY-CKAB*G1-(PI/2)*BEFF*CKAB*G1)

C	   CKAB=(FY-Q)/G1

	   ZETAB=BEFF
        CDAB=2.0*TM*ZETAB*WB

        ALP=FY/(G1*CKAB)
        ALP=1.0/ALP
       
        
C       ALP=(FY-Q)/FY
	  
        DT=DTN/NDIV
        WRITE(3,1321)T1
        WRITE(3,1336)RMBM        
        WRITE(3,1301)TB
        WRITE(3,1337)BM
        WRITE(3,1304)ZETAB
        WRITE(3,1342)CKAB
        WRITE(3,1345)BEFF
        WRITE(3,1344)FY
        WRITE(3,1339)NDIV
C       COMPUTATION OF SUPERSTRUCTURE MATRICES  
        CALL SMS(SM,AM,NST)
        CALL SKS(SK,VI,AK,NST)
        CALL INVD(VI,SM,V,NST)
        CALL EIGEN(VI,SM,V,W,NST,99)
        IF(KNOR.EQ.0) GO TO 102
        CALL NORM(W(1),W1,AK,NST)
        CALL SKS(SK,VI,AK,NST)
	   WRITE(*,*)((SK(I,J),J=1,NST),I=1,NST)

        CALL INVD(VI,SM,V,NST)
        CALL EIGEN(VI,SM,V,W,NST,99)
102     CALL INV(V,VI,NST)
        CALL SCS(CD,VI,W,ZETA,NST)
        WRITE(3,1351)
        WRITE(3,1325)
        WRITE(3,1302)((SM(I,J),J=1,NST),I=1,NST)
        WRITE(3,1300)
        WRITE(3,1302) (W(I),I=1,NST)
        WRITE(3,1329)
        WRITE(3,1302)((V(I,J),J=1,NST),I=1,NST)
        WRITE(3,1326)
        WRITE(3,1302)((CD(I,J),J=1,NST),I=1,NST)
        WRITE(3,1310)
        WRITE(3,1302)((SK(I,J),J=1,NST),I=1,NST)
        CALL SMO(SM,BM,NDOF)
        CALL SCO(CD,CDAB,C1,NDOF)
        CALL SKO(SK,CKAB,AK1,NDOF)
        WRITE(3,1328)
        WRITE(3,1312)((SM(I,J),J=1,NDOF),I=1,NDOF)
        WRITE(3,1327)
        WRITE(3,1312)((CD(I,J),J=1,NDOF),I=1,NDOF)
        WRITE(3,1311)
        WRITE(3,1312)((SK(I,J),J=1,NDOF),I=1,NDOF)
C       INITIALIZATION OF VECTORS.
        DO 1030 KW=1,NW
        DO 110 I=1,NDOF
        D1(I)=0.0
        V1(I)=0.0
        A1(I)=0.0
        D2(I)=0.0
        V2(I)=0.0
        A2(I)=0.0
        AA(I)=0.0
        PD(I)=0.0
        PA(I)=0.0
        DD(I)=0.0
        DV(I)=0.0
        P1(I)=0.0
        P2(I)=0.0
        PINC(I)=0.0
 110    CONTINUE
        TIME=0.0
        FAB=0.0
        PZ=0.0
        Z=0.0
        WRITE(8,*) 'I AM HERE'
        DO 1020 K=1,NDT-1
        IF(KSN.EQ.1) XG(K)=AMP*SIN(WW(KW)*TIME)
        TT=TIME+DTN
        IF(KSN.EQ.1) XG(K+1)=AMP*SIN(WW(KW)*TT)
        DO 112 I=1,NDOF
        P1(I)=-SM(I,I)*XG(K)
        P2(I)=-SM(I,I)*XG(K+1)
 112    CONTINUE
        DO 116 I=1,NDOF
        PINC(I)=(P2(I)-P1(I))/NDIV
 116    CONTINUE
        DO 1010 KS=1,NDIV
        TIME = TIME + DT
        DO 117 I=1,NDOF
        P2(I)=P1(I)+PINC(I)
 117    CONTINUE        
        EQ=-P2(1)/SM(1,1)
        IF(K+KS-2)1020,118,120       
 118    WRITE(8,*) (P2(I),I=1,NDOF)
        WRITE(8,*) (P1(I),I=1,NDOF)
        write(8,*) (PINC(I),I=1,NDOF)
        WRITE(8,*) PZ
        CALL ACCN(P2,NDOF,PZ) 
        DO 119 I=1,NST
        A1(I)=A2(I)
        A2(I)=0.0
 119    CONTINUE
        WRITE(8,*) (A1(I),I=1,NDOF)
 120    DPZ=0.0
        DO 1005 KIT=1,NIT
        WRITE(8,*)DPZ
        CALL EFFP(PINC,EFP,DT,NDOF,DPZ)
        WRITE(8,*) (EFP(I),I=1,NDOF)
        CALL EFFK(EFK,DT,NDOF)
        WRITE(8,*) ((EFK(I,J),J=1,NDOF),I=1,NDOF)
        CALL GAUS(NDOF,EFK,EFP,DD)
        WRITE(8,*)(DD(I),I=1,NDOF)
        
        DO 142 I=1,NDOF
        D2(I)=D1(I)+DD(I)
        DV(I)=3.0*DD(I)/DT-3.0*V1(I)-DT*A1(I)*0.5
        V2(I)=V1(I)+DV(I)
 142    CONTINUE
 143    CALL WEN(V2(NDOF),Z,DZ,DPZ)
        write(8,*) (K-1)*NDIV+KS, KIT, DPZ, DZ  
        WRITE(8,*)'o-o-o-o-o-o'
 1005   CONTINUE
        PZ=PZ+DPZ
        Z=Z+DZ
        FAB=CKAB*D2(NDOF)+PZ
        FAB1 = 0.0
        
        CALL ACCN(P2,NDOF,PZ)
        CALL PEAK(P1,P2,AA,EQ,NDOF)
 1010   CONTINUE
        RF=0.0
        DO 151 I=1,NDOF
        RF=RF - SM(I,I)*AA(I)
 151    CONTINUE
        WRITE(7,*) 'Break'
        DO 2424 KK = 1,NST
        FAB1 = FAB1 - SM(KK,KK)*AA(KK)
        WRITE(7,*) SM(KK,KK)
        WRITE(7,*) AA(KK)
2424    CONTINUE
C        WRITE(*,1316) K,TIME,FAB,RF
        WRITE(4,1346) TIME, EQ, AA(1),V2(1),D2(NDOF), FAB, FAB1
C	1  ,FAB
 1020   CONTINUE
        
C	  WRITE(*,*)ABS(PD(NDOF))*100

	  AQ=Q/(TM*9.81)
	  BKB=(2*CKAB*ABS(PD(NDOF)))/(TM*9.81)
	  CA=TM*ABS(PA(1))/(TM*9.81)
	  RE=AQ+BKB+CA

	  WRITE(5,1356) IJK,ABS(PA(1)),ABS(PD(1))
	 
       WRITE(3,1350)
        DO 152 I=1,NDOF
        WRITE(3,1349) PD(I),PA(I)
 152    CONTINUE
 1030   CONTINUE
 1050   CONTINUE
C
 1300   FORMAT(/,10X,'FIXED BASE NATURAL FREQUENCIES ')
 1301   FORMAT(3X,'PERIOD OF ISOLATED STRUCTURE     = ',E13.5)
 1302   FORMAT(5X,5E13.6)
 1303   FORMAT(/10X,'STIFFNESS MARTRIX')
 1304   FORMAT(3X,'DAMPING RATIO OF SLIDING SYSTEM  = ',E13.5)
 1305   FORMAT(/10X,'NATURAL FREQUENCIES')
 1306   FORMAT(/10X,'MODAL COLUMN MARTRIX')
 1310   FORMAT(/,10X,'FIXED BASE STIFFNESS MATRIX')
 1311   FORMAT(/,10X,'ASSEMBLED STIFFNESS MATRIX')
 1312   FORMAT(5X,6E13.5)
 1313   FORMAT(3X,4E13.5)
 1315   FORMAT(6X,'NDT   TIME       FABX  ',8X,'FABY')
 1316   FORMAT(5X,I10,F8.3,5X,3E15.6)
 1319   FORMAT(3X,I2,2F9.2,2X,2F10.2)
 1321   FORMAT(3X,'FIRST SUPERSTRUCTURE PERIOD      = ',E13.5)
 1325   FORMAT(/,10X,'FIXED BASE MASS MATRIX')
 1326   FORMAT(/,10X,'FIXED BASE DAMPING MATRIX')
 1327   FORMAT(/,10X,'ASSEMBLED DAMPING MATRIX')
 1328   FORMAT(/,10X,'ASSEMBLED MASS MATRIX')
 1329   FORMAT(/,10X,'FIXED BASE MODE-SHAPE MATRIX')
 1336   FORMAT(3X,'RATION OF BASE MASS TO FLOOR     = ',E13.5)
 1337   FORMAT(3X,'MASS OF BASE RAFT                = ',E13.5)
 1339   FORMAT(3X,'NDIV VALUE                       = ',I13)
 1342   FORMAT(3X,'SLIDING SYSTEM STIFFNESS         = ',E13.5)
 1344   FORMAT(3X,'YIELD FORCE LEVEL                = ',E13.5)
 1345   FORMAT(3X,'FY BY W RATIO                    = ',E13.5)
 1346   FORMAT(F6.3,1X,7E13.5)
 1349   FORMAT(3X,4E13.5)
 1350   FORMAT(/,5X,'PEAK DISP    PEAK ACCN      ( X/Y)')
 1351   FORMAT(/,3X,'SUPERSTRUCTURE PROPERTIES  =  ',15(' X'))
 1353   FORMAT(I3,6X,6F6.2,F6.2,F6.2)
 1354   FORMAT(6F10.4)
 1356   FORMAT(I10,6E13.5)
 1362   FORMAT(3X,6F7.2,4F8.2,2F5.2)
 1367   FORMAT(I4,F10.2,6E14.6)
        STOP
        END
 
        SUBROUTINE SMS(SM,AM,NST)
        DIMENSION SM(11,11),AM(11)
        DO 1 I=1,NST+1
        DO 1 J=1,NST+1
        SM(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,NST
        SM(I,I)=SM(I,I)+AM(I)
   2    CONTINUE
        RETURN
        END
 
        SUBROUTINE SKS(SK,VI,AK,NST)
        DIMENSION SK(11,11),VI(11,11),AK(11)
        DO 1 I=1,NST+1
        DO 1 J=1,NST+1
        SK(I,J)=0.0
   1    CONTINUE
        SK(1,1)=SK(1,1) + AK(1)
        SK(1,2)=SK(1,2) - AK(1)
        SK(2,1)=SK(2,1) - AK(1)
        DO 2 I=2,NST
        II=I-1
        SK(I,I)=AK(I)+AK(II)
        IF(I.EQ.NST) GO TO 2
        SK(I,I+1)=SK(I,I+1)-AK(I)
        SK(I+1,I)=SK(I+1,I)-AK(I)
   2    CONTINUE
        DO 3 I=1,NST
        DO 3 J=1,NST
        VI(I,J)=SK(I,J)
   3    CONTINUE
        RETURN
        END
 
        SUBROUTINE SCS(CD,VI,W,ZETA,NST)
        DIMENSION CD(11,11),VI(11,11),W(11),ZETA(11),XX(11,11)
        DO 1 I=1,NST+1
        DO 1 J=1,NST+1
        CD(I,J)=0.0
   1    CONTINUE
        DO 2 I=1,NST
        CD(I,I)=2.0*ZETA(I)*W(I)
   2    CONTINUE
        DO 3 I=1,NST
        DO 3 J=1,NST
        XX(I,J)=0.0
        DO 3 K=1,NST
        XX(I,J)=XX(I,J) + CD(I,K)*VI(K,J)
   3    CONTINUE
        DO 4 I=1,NST
        DO 4 J=1,NST
        CD(I,J)=0.0
        DO 4 K=1,NST
        CD(I,J)=CD(I,J) + VI(K,I)*XX(K,J)
   4    CONTINUE
        RETURN
        END
 
        SUBROUTINE SMO(SM,BM,NDOF)
        DIMENSION SM(11,11)
        NST=NDOF-1
        DO 1 I=1,NST
        SM(I,NST+1)=SM(I,I)
   1    CONTINUE
        SM(NDOF,NDOF)=BM
        RETURN
        END
 
        SUBROUTINE SCO(CD,CDAB,C1,NDOF)
        DIMENSION CD(11,11)
        NST=NDOF-1
        C1=CD(NST,NST)+CD(NST,NST-1)
        CD(NDOF,NDOF)=CDAB
        CD(NDOF,NST)=-C1
        RETURN
        END
 
        SUBROUTINE SKO(SK,CKAB,AK1,NDOF)
        DIMENSION SK(11,11)
        NST=NDOF-1
        AK1=SK(NST,NST)+SK(NST,NST-1)
        SK(NDOF,NDOF)=CKAB
        SK(NST,NDOF)=0.0
        SK(NDOF,NST)=-AK1
        RETURN
        END
 
        SUBROUTINE EFFP(PINC,EFP,T,NDOF,DPZ)
        COMMON/RS5/SM,SK,CD
        COMMON/RS6/D1,V1,A1
        DIMENSION EFP(11),PINC(11),V1(11),A1(11)
        DIMENSION SM(11,11),SK(11,11),CD(11,11),D1(11)
        DO 1 I=1,NDOF
   1    EFP(I) =0.0
        DO 2 I=1,NDOF
        EFP(I)=EFP(I)+PINC(I)
        DO 3 J=1,NDOF
        X1=SM(I,J)*(6.0*V1(J)/T + 3.0*A1(J))
        X2=CD(I,J)*(3.0*V1(J) + 0.5*T*A1(J))
        EFP(I)=EFP(I) + X1 + X2
   3    CONTINUE
   2    CONTINUE
        EFP(NDOF)=EFP(NDOF)-DPZ
C       WRITE(8,11)
C       WRITE(8,12)(EFP(I),I=1,3)
C 11    FORMAT(10X,'EFFECTIVE FORCE')
        RETURN
        END
 
        SUBROUTINE EFFK(EFK,T,NDOF)
        DIMENSION EFK(11,11)
        DIMENSION SK(11,11),SM(11,11),CD(11,11)
        COMMON/RS5/SM,SK,CD
        DO 1 I=1,NDOF
        DO 1 J=1,NDOF
   1    EFK(I,J)=SK(I,J)
        DO 2 I=1,NDOF
        DO 2 J=1,NDOF
        EFK(I,J)=EFK(I,J) + CD(I,J)*3.0/T + SM(I,J)*6.0/(T*T)
   2    CONTINUE
C       WRITE(8,11)
C       WRITE(8,12)((EFK(I,J),J=1,NDOF),I=1,NDOF)
C 11    FORMAT(10X,'EFFECTIVE STIFFNESS MATRIX')
C 12    FORMAT(5X,3E13.6)
        RETURN
        END
 
        SUBROUTINE GAUS(N,A,F,X)
        DIMENSION A(11,11),F(11),X(11)
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
 30     CONTINUE
 20     CONTINUE
 10     CONTINUE
        X(N)=F(N)/A(N,N)
        DO 40 I=(N-1),1,-1
        SUM=0.0
        DO 50 J=(I+1),N
        SUM=SUM+A(I,J)*X(J)
 50     CONTINUE
        X(I)=(F(I)-SUM)/A(I,I)
 40     CONTINUE
C       WRITE(*,11)
C       WRITE(*,12)(X(I),I=1,NDOF)
C 11    FORMAT(10X,'DELTA DISPLACEMENT')
C 12    FORMAT(3X,2E13.6)
        RETURN
        END

        SUBROUTINE ACCN(P,NDOF,F)
        DIMENSION SK(11,11),SM(11,11),CD(11,11)
        DIMENSION FK(11),FC(11),EP1(11),P(11)
        DIMENSION D(11),V(11),A(11)
        COMMON/RS5/SM,SK,CD
        COMMON/RS7/D,V,A
        DO 3 I=1,NDOF
        FC(I)=0.0
        FK(I)=0.0
        DO 2 J=1,NDOF
        FC(I)=FC(I)+CD(I,J)*V(J)
        FK(I)=FK(I)+SK(I,J)*D(J)
  2     CONTINUE
        EP1(I)=P(I)-FC(I)-FK(I)
  3     CONTINUE
        EP1(NDOF)=EP1(NDOF)-F
        A(NDOF)=EP1(NDOF)/SM(NDOF,NDOF)
        DO 4 I=1,NDOF-1
        A(I)=(EP1(I) - SM(I,I)*A(NDOF))/SM(I,I)
  4     CONTINUE
        RETURN
        END

        SUBROUTINE PEAK(P1,P2,AA,EQ,N)
        DIMENSION D1(11),A1(11),PD(11),V1(11),PV(11)
        DIMENSION D2(11),A2(11),PA(11),V2(11)
        DIMENSION P1(11),P2(11),AA(11) 
        COMMON/RS6/D1,V1,A1 
        COMMON/RS7/D2,V2,A2
        COMMON/RS8/PD,PV,PA
        AA(N)=A2(N)+EQ
        DO 1 I=1,N-1         
        AA(I)=A2(I)+A2(N)+EQ
   1    CONTINUE
        DO 8 I=1,N
        IF(ABS(PD(I))-ABS(D2(I)))7,8,8
   7    PD(I)=D2(I)
   8    CONTINUE
        DO 10 I=1,N
        IF(ABS(PV(I))-ABS(V2(I)))9,10,10
   9    PV(I)=V2(I)
  10    CONTINUE
        DO 14 I=1,N
        IF(ABS(PA(I))-ABS(AA(I)))13,14,14
  13    PA(I)=AA(I)
  14    CONTINUE
        DO 17 I=1,N
        A1(I)=A2(I)
        V1(I)=V2(I)
        D1(I)=D2(I)
        P1(I)=P2(I)
  17    CONTINUE
        RETURN
        END

        SUBROUTINE EIGEN(XK,XM,V,W,N,IT)
C       STODOLA METHOD FOR MODE SHAPE AND FREQUENCY.
        DIMENSION XK(11,11),S(11,11),X(100,11),X1(11),X2(11)
        DIMENSION XX(11,11),XM(11,11),V(11,11),W(11)
        DO 1 I=1,N
        DO 1 J=1,N
        S(I,J)=0.0
        S(I,I)=1.0
   1    CONTINUE
        DO 16 IE=1,N
        IF(IE-1) 16,6,2
   2    DO 4 I=1,N
        DO 4 J=1,N
        XX(I,J)= 0.0
        DO 3 K=1,N
   3    XX(I,J)=XX(I,J)+X1(I)*X1(K)*XM(K,J)
   4    S(I,J)=S(I,J)-XX(I,J)/UM
        DO 5 I=1,N
        DO 5 J=1,N
        XK(I,J)=0.0
        DO 5 K=1,N
   5    XK(I,J)=XK(I,J)+X(I,K)*S(K,J)
   6    DO 7 I=1,N
   7    X(1,I)=1.0
        DO 11 K=1,IT
        DO 9 I=1,N
        X2(I)=0.0
        DO 8 J=1,N
   8    X2(I)=X2(I)+XK(I,J)*X(K,J)
   9    X(K+1,I)=X2(I)/X2(1)
  11    CONTINUE
        W(IE)=SQRT(1.0/X2(1))
        DO 12 I=1,N
        X1(I)=X(IT,I)
        DO 12 J=1,N
  12    X(I,J)=XK(I,J)
        UM=0.0
        DO 14 I=1,N
        X2(I) = 0.0
        DO 13 J = 1,N
  13    X2(I)=X2(I)+XM(I,J)*X1(J)
  14    UM=UM+X1(I)*X2(I)
        DO 15 I=1,N
  15    V(I,IE)=X1(I)/SQRT(UM)
  16    CONTINUE
        RETURN
        END
 
        SUBROUTINE INVD(XK,XM,V,N)
C       GENERATION OF FLEXIBILITY & DYNAMIC MATRIX.
        DIMENSION XK(11,11),V(11,11),XM(11,11)
        DO 1 I=1,N
        DO 1 J=1,N
        V(I,J)=0.0
        V(I,I)=1.0
   1    CONTINUE
        DO 6 I=1,N
        TEMP=XK(I,I)
        DO 2 J=1,N
        XK(I,J)=XK(I,J)/TEMP
        V(I,J)=V(I,J)/TEMP
   2    CONTINUE
        DO 5 J=1,N
        IF(I-J) 3,5,3
   3    RT=XK(J,I)
        DO 4 K=1,N
        XK(J,K)=XK(J,K)-XK(I,K)*RT
   4    V(J,K)=V(J,K)-V(I,K)*RT
   5    CONTINUE
   6    CONTINUE
        DO 7 I=1,N
        DO 7 J=1,N
        XK(I,J)=0.0
        DO 7 K=1,N
   7    XK(I,J)=XK(I,J)+V(I,K)*XM(K,J)
        RETURN
        END
		 
        SUBROUTINE INV(XX,B,N)
C       TO INVERT A REAL MATRIX.
        DIMENSION A(11,11),B(11,11),XX(11,11)
        DO 1 I=1,N
        DO 1 J=1,N
        A(I,J)=XX(I,J)
        B(I,J)=0.0
        B(I,I)=1.0
   1    CONTINUE
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

        SUBROUTINE NORM(W1,WS1,AK,NST)
        DIMENSION AK(11)
        RT=WS1/W1
        DO 11 I=1,NST
        AK(I)=AK(I)*RT*RT
  11    CONTINUE
        RETURN
        END

        SUBROUTINE WEN(V,Z,DZ,DPZ)
        COMMON/RS12/FY,G1,ALP,BT,G2,A,NT,DT

        B1=DT*(-G2*ABS(V)*Z*(ABS(Z)**(NT-1)))
     .    +DT*(-BT*V*(ABS(Z)**NT) + A*V)

        B2=DT*(-G2*ABS(V)*(Z+B1*0.5)*(ABS(Z+B1*0.5)**(NT-1)))
     .  +DT*(-BT*V*(ABS(Z+B1*0.5)**NT) + A*V)

        B3=DT*(-G2*ABS(V)*(Z+B2*0.5)*(ABS(Z+B2*0.5)**(NT-1)))
     .    +DT*(-BT*V*ABS((Z+B2*0.5)**NT) + A*V)

        B4=DT*(-G2*ABS(V)*(Z+B3)*(ABS(Z+B3)**(NT-1)))
     .    +DT*(-BT*V*(ABS(Z+B3)**NT) + A*V)

        DZ=(B1+(B2+B3)*2.0+B4)/(6.0*G1)
        
        DPZ=(1.0-ALP)*FY*DZ
        
        RETURN
        END
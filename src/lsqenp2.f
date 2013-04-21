      SUBROUTINE LSQENP(N,K,M,Y,X,B,IP,IB,IDVT,ICON,IQUIT,IPRNT)
C LSQENP REVISED JUNE 24, 1976
C LSQENP AND CONLIM REVISED FEB., 1975
C REVISED 9-18-73
C  REVISED 11-20-72 TO ADD OPTION FOR NOT PRINTING OBSERVED AND CALCULAT
C      VALUES AFTER CONVERSION.  SET IPRNT LESS THAN 0.
C     REVISED 6-12-70 TO CALCULATE SUM(( F -Y(I))**2) IN DBLE PRECISION
C        NONLINEAR LEAST SQUARES
      CHARACTER*1 JBCH,JOCH,JPCH,JXCH,JYCH
      DIMENSION BS(50),DB(50),BA(50),G(50),SA(50),P(50),A(50,51)
      DIMENSION X(*),Y(*),B(*),IB(*)
      COMMON/LSQPLT/IFP,YMIN,YMAX
      COMMON /PARAM/FF,T,E,TAU,AL,GAMCR,DEL,ZETA
      COMMON/WEITS/WTS(500000)
      COMMON /AFCLSQ/ SE,PHI,PHIZ,WS,XL,IFSS2,IFSS3,IWS6,I,J,JJ,IBK2,
     1                IBKA,IBKM,IPLOT,BS,SA,A
      REAL*8 PHD,XLL,DTG,GTG    

C     Threadprivate directive to make common blocks local 
C     and persistent to a thread
#ifdef _OPENMP
!$OMP THREADPRIVATE(/LSQPLT/,/PARAM/,/WEITS/,/AFCLSQ/)
#endif
C     open(26,file='o_lsqenp')
      JBCH=' '
      JOCH='O'
      JPCH='P'
      JXCH='X'
      JYCH='Y'
      IWS4=IQUIT
      IWS6=ICON
      IPLOT=IFP
      IF (IFP .EQ. 1) SPRD=YMAX-YMIN
C     ------------------------------------------------------------------
C     MAX NO OF PARAMETERS IS K=50
C     MAX NO OF OBSERVATIONS IS N=500000
C                  IWHER = 1 MEANS GET P S AND F
C                  IWHER GREATER THAN 1 MEANS GET F ONLY
   23 IWHER = 0
      GO TO 43
   31 CONTINUE
      CALL FCODE(Y,X,B,F,I,N)
      IF (IWHER.NE.1)GO TO 39
   34 IF (IFSS2.NE.0)GO TO 39
   35 CONTINUE
      CALL PCODE(P,X,B,F,I)
   39 CONTINUE
      IF (IWHER.EQ.0)GO TO 43
C             1  2   3   4
   42 GO TO(106,363,114,125), IWHER
   43 ITCT=0
      IF (IPLOT .LE. 0) GO TO 55
   49 CONTINUE
c      IBCH=JBCH
c      IOCH=JOCH
c      IPCH=JPCH
c      IXCH=JXCH
c      IYCH=JYCH
      IBCH=ichar(JBCH)
      IOCH=ichar(JOCH)
      IPCH=ichar(JPCH)
      IXCH=ichar(JXCH)
      IYCH=ichar(JYCH)
   55 IF(IP.LE.0)GOTO 62
      DO 61 I=1,IP
      IF (IB(I).GT.0)GO TO 61
C  58 WRITE (26,402)
   58 CONTINUE
      STOP
   61 CONTINUE
   62 CONTINUE
      XKDB = 1.
C     ..................................................................
C                  START THE CALCULATION OF THE PTP MATRIX
      IBKA=1
      IF(IPRNT .LE. 0) GO TO 82
C     WRITE (26,383)N,K,IP,M,IFP,GAMCR,DEL,FF,T,E,TAU,AL,ZETA
   82 CONTINUE
      DO 86 I=1,K
      G (I) =0.
      DO 86 J=1,K
   86 A (I,J)=0.
      GO TO(88,91,91),IBKA
   88 IFSS3=IPRNT
      IFSS2=IDVT
      GO TO 92
   91  IFSS3=1
C  92 IF(IFSS3 .GT. 0) WRITE(26,384)(B(J),J=1,K)
   92 CONTINUE
C  92 WRITE(26,384)(B(J),J=1,K)
       IF (IFSS3 .LE. 0 .OR. IPRNT.LT. 0) GO TO 99
      IF (IPLOT .LE. 0) GO TO 98
      WS =YMIN+SPRD
C     WRITE(26,382)YMIN,WS
      GO TO 99
C  98 WRITE (26,386)
   98 CONTINUE
   99 I=1
      PHD=0.
      IF (IFSS2.EQ.0)GO TO 104
      GO TO 111
  103 IF (IFSS2.EQ.1)GO TO 112
  104 IWHER=1
C                  GET P S AND F
      GO TO 31
  106 IF (IP.LE.0)GO TO 132
  107 DO 109 II=1,IP
      IWS=IB(II)
  109 P(IWS)=0.
      GO TO 132
C     ..................................................................
C                  THIS IS THE ESTIMATED P S ROUTINE
  111 CONTINUE
  112 IWHER=3
      GO TO 31
  114 FWS=F
      J=1
  116 IF (IP.LE.0)GO TO 120
  117 DO 119 II=1,IP
      IF ((J-IB(II)).EQ.0)GO TO 128
  119 CONTINUE
  120 DBW=B(J)*DEL
       IF (B(J) .EQ. 0.) DBW = DEL
      TWS=B(J)
      B(J)=B(J)+DBW
      IWHER=4
      GO TO 31
  125 B(J)=TWS
      P(J)=(F-FWS)/DBW
      GO TO 129
  128 P(J)=0.
  129 J=J+1
      IF ((J-K).LE.0)GO TO 116
  131 F=FWS
C                  END OF ESTIMATED P S ROUTINE
C     ..................................................................
C                  NOW, USE THE P S TO MAKE PARTIALS MATRIX
  132 DO 136 JJ=1,K
      G(JJ)=G(JJ)+(Y(I)-F)*P(JJ)
      DO 136 II = JJ,K
      A(II,JJ)=A(II,JJ)+P(II)*P(JJ)
  136 A(JJ,II)=A(II,JJ)
      IF (IPLOT .LE. 0) GO TO 184
  138 IF (IFSS3.LE.0)GO TO 188
C                  PLOTTING Y(I),F
  139 IO = (Y(I)-YMIN)*100./SPRD
      IPP = (F-YMIN)*100./SPRD
      IF (IO.EQ.IPP)GO TO 148
      IF (IO.GT. IPP)GO TO 153
C                  Y(I) OUT FIRST
  143 IP1=IOCH
      IP2=IPCH
      I1=IO
      I2=IPP
      GO TO 157
C                  ONLY ONE CHARACTER
  148 IP1=IYCH
      IP2=IBCH
      I1=IO
      I2=IPP
      GO TO 157
C                  F OUT FIRST
  153 IP1=IPCH
      IP2=IOCH
      I1=IPP
      I2=IO
C                  ZERO PLOTS IN THE LEFT HAND COLUMN, SO I1 IS ITS
C                  OWN BLANK COUNTER
C                  OVERFLOWS PLOT X IN COLUMN 102
C                  UNDERFLOWS ALSO PLOT X IN COLUMN ZERO
  157 IF (I2.LE.101)GO TO 165
  158 I2=101
      IP2=IXCH
      IF (I1.LT.101)GO TO 165
  161 I1=101
      IP1=IXCH
      IP2=IBCH
      GO TO 171
  165 IF (I1.GE.0)GO TO 171
  166 I1=0
      IP1=IXCH
      IF (I2.GT.0)GO TO 171
  169 I2=1
      IP2=IBCH
  171 I1M1=I1
      I1M2=I2-I1-1
      IF (I1M1.GT.0)GO TO 179
  174 IF (I1M2.GT.0)GO TO 177
C 175 WRITE (26,404)IP1,IP2
  175 CONTINUE
      GO TO 188
C 177 WRITE (26,404)IP1,(IBCH,II=1,I1M2),IP2
  177 CONTINUE
      GO TO 188
  179 IF (I1M2.GT.0)GO TO 182
C 180 WRITE (26,404)(IBCH,II=1,I1M1),IP1,IP2
  180 CONTINUE
      GO TO 188
C 182 WRITE (26,404)(IBCH,II=1,I1M1),IP1,(IBCH,II=1,I1M2),IP2
  182 CONTINUE
      GO TO 188
  184 WS=Y(I)-F
       IF (IFSS3 .LE. 0 .OR. IPRNT.LT. 0) GO TO 188
C 187  WRITE (26,401)X(I),Y(I),F,WS      
  187  CONTINUE
  188 WS=WTS(I)*(Y(I)-F)
      PHD=PHD+WS*WS*1.0D0
      I=I+1
      IF (I.LE.N)GO TO 103
      PHI=PHD
      IF (IP.LE.0)GO TO 199
  193 DO 198 JJ=1,IP
      IWS=IB(JJ)
      DO 197 II=1,K
      A(IWS,II)=0.
  197 A(II,IWS)=0.
  198 A(IWS,IWS)=1.
C     INSERT 1
  199 IF(IBKA .EQ. 1)GO TO 204
      IBKS=1
      CALL       CONLIM(N,K,M,Y,X,B,IP,IB,IBKS,IBD)
      GO TO(82,314,38,359,38),IBD
C                  SAVE SQUARE ROOTS OF DIAGONAL ELEMENTS
  204 CONTINUE
      DO 206 I=1,K
  206 SA(I)=SQRT (A(I,I))
      DO 219 I=1,K
      DO 214 J=1,K
      WS = SA(I)*SA(J)
      IF(WS.GT.0.)GOTO 213
      A(I,J) =0.
      GO TO 214
  213 A(I,J)=A(I,J)/WS
  214 CONTINUE
      IF(SA(I).GT.0.)GOTO 218
      G(I)=0.
      GO TO 219
  218 G(I)=G(I)/SA(I)
  219 CONTINUE
      DO 221 I=1,K
  221 A(I,I)=1.
      PHIZ=PHI
C                  WE NOW HAVE PHI ZERO
C     ..................................................................
      IF (ITCT.GT.0)GO TO 230
C                  FIRST ITERATION
  225 XL=AL
      ITCT=1
      DO 229 J=1,K
  229 BS(J)=B(J)
C                  BS(J) CORRESPONDS TO PHIZ
  230 IBK1=1
      WS=N-K+IP
      SE=SQRT(PHIZ/WS)
      IF (IFSS3.GT.0)GO TO 239
      IF(IPRNT .LE. 0) GO TO 310
  234 IF (IFSS2.EQ.0)GO TO 237
C 235 WRITE (26,387)PHIZ,SE,XLL,GAMMA,XL
  235 CONTINUE
      GO TO 310
C 237 WRITE (26,388)PHIZ,SE,XLL,GAMMA,XL
  237 CONTINUE
      GO TO 310
  239 IF (IFSS2.EQ.0)GO TO 242
C 240 WRITE (26,379)PHIZ,SE,XL
  240 CONTINUE
      GO TO 310
C 242 WRITE (26,385)PHIZ,SE,XL
  242 CONTINUE
      GO TO 310
  244 PHIL=PHI
C                  WE NOW HAVE PHI LAMBDA
      DO 247 J=1,K
      IF (ABS(DB(J)/(ABS(B(J)) + TAU)).GE.E)GOTO 251
  247 CONTINUE
C     WRITE (26,399)
      IBS=4
      GO TO 371
  251 IF (IWS4.EQ.0)GO TO 257
      IF (IWS4.EQ.1)GO TO 255
      IWS4=IWS4-1
      GO TO 257
C 255 WRITE (26,400)
  255 CONTINUE
      GO TO 371
  257 XKDB = 1.
      IF (PHIL.GT.PHIZ)GO TO 281
  259 XLS=XL
      DO 262 J=1,K
      BA(J)=B(J)
  262 B(J)=BS(J)
      IF (XL.GT..00000001)GO TO 268
  264 DO 266 J=1,K
      B(J)=BA(J)
  266 BS(J)=B(J)
      GO TO 82
  268 XL=XL/10.
      IBK1=2
      GO TO 310
  271 PHL4=PHI
C                  WE NOW HAVE PHI(LAMBDA/10)
      IF(PHL4.GT.PHIZ)GOTO 276
  273 DO 274 J=1,K
  274 BS(J)=B(J)
      GO TO 82
  276 XL=XLS
      DO 279 J=1,K
      BS(J)=BA(J)
  279 B(J)=BA(J)
      GO TO 82
  281 IBK1=4
      XLS=XL
      XL=XL/10.
      DO 285 J=1,K
  285 B(J)=BS(J)
      GO TO 310
  287 IF (PHI.LE.PHIZ)GO TO 296
  288 XL=XLS
      IBK1=3
  290 XL=XL*10.
  291 DO 292 J=1,K
  292 B(J)=BS(J)
      GO TO 310
  294 PHIT4=PHI
C                  WE NOW HAVE PHI(10*LAMBDA)
      IF (PHIT4.GT.PHIZ)GO TO 299
  296 DO 297 J=1,K
  297 BS(J)=B(J)
      GO TO 82
  299 IF (GAMMA.GE.GAMCR)GO TO 290
  300 XKDB = XKDB/2.
      DO 303 J=1,K
      IF (ABS(DB(J)/(ABS(B(J))+TAU)).GE.E)GO TO 291
  303 CONTINUE
      DO 305 J=1,K
  305 B(J)=BS(J)
C     WRITE (26,410)
      IBS=4
      GO TO 371
C
C     ..................................................................
C                  SET UP FOR MATRIX INVERSION
  310 CONTINUE
      DO 312 I=1,K
  312 A(I,I)=A(I,I)+XL
C                  GET INVERSE OF A AND SOLVE FOR DB (J)S
      IBKM=1
C     ..................................................................
C                  THIS IS THE MATRIX INVERSION ROUTINE
C                  K IS THE SIZE OF THE MATRIX
  314 CALL GJR(A,K,ZETA,MSING)
      GO TO(316,38),MSING
C     INSERT 2
  316 IF (IBKM .EQ. 1)GO TO 321
      IBKS=2
      CALL       CONLIM(N,K,M,Y,X,B,IP,IB,IBKS,IBD)
      GO TO(82,314,38,359,38),IBD
C                  END OF MATRIX INVERSION, SOLVE FOR DB(J)
  321 DO 325 I=1,K
      DB(I)=0.
      DO 324 J=1,K
  324 DB(I)=A(I,J)*G(J)+DB(I)
  325 DB(I)=XKDB*DB(I)
      XLL=0.0
      DTG = 0.
      GTG = 0.
      DO 334 J=1,K
      XLL=XLL+DB(J)*DB(J)
      DTG = DTG + DB(J)*G(J)
      GTG = GTG + G(J)**2
      IF (SA(J) .GT. 0.0) GO TO 333
C     WRITE (26,332) J,SA(J)
      GO TO 335
  333 CONTINUE
      DB(J)=DB(J)/SA(J)
  335 CONTINUE
      B(J)=B(J)+DB(J)
  334 CONTINUE
      KIP=K-IP
      IF (KIP.EQ.1)GO TO 350
      CGAM=DTG/DSQRT(XLL*GTG)
      JGAM = 1
      IF(CGAM.GT..0)GOTO 342
  340 CGAM = ABS(CGAM)
      JGAM = 2
  342 IF(CGAM .GT. 1.0)CGAM=1.0
      GAMMA = 57.2957795*(1.5707288+CGAM*(-0.2121144+CGAM*(0.074261
     1-CGAM*.0187293)))*SQRT(1.-CGAM)
      GO TO(351,344), JGAM
  344 GAMMA = 180.-GAMMA
      IF (XL.LT.1.0)GO TO 351
C 346 WRITE(26,398)XL,GAMMA
  346 CONTINUE
      IBS=4
      GO TO 371
  350 GAMMA=0.
  351 XLL=DSQRT(XLL)
      IBK2=1
      GO TO 359
  354 IF (IFSS3.LE.0)GO TO 358
C 355 WRITE (26,380)(DB(J),J=1,K)
  355 CONTINUE
C     WRITE (26,381)PHI,XL,GAMMA,XLL
  358 GO TO(244,271,294,287),IBK1
C
C     ..................................................................
C                  CALCULATE PHI
  359 I=1
      PHD=0.
      IWHER=2
      GO TO 31
  363 CONTINUE
      PHD=PHD+((Y(I)-F)*1.0D0*WTS(I))**2
      I=I+1
      IF (I.LE.N)GO TO 31
      PHI=PHD
C     INSERT 3
      IF (IBK2 .EQ. 1)GO TO 354
      IBKS=3
      CALL       CONLIM(N,K,M,Y,X,B,IP,IB,IBKS,IBD)
      GO TO(82,314,38,359,38),IBD
  371 CONTINUE
      IBKS=4
      CALL       CONLIM(N,K,M,Y,X,B,IP,IB,IBKS,IBD)
      GO TO(82,314,38,359,38),IBD
C     ..................................................................
   38 RETURN
  332 FORMAT(5X,'J=',I5,5X,'SA(J) = ',E15.6)
  376 FORMAT (25I3)
  377 FORMAT (7F10.0)
  378 FORMAT(12A6)
  379 FORMAT (/,13X,' PHI',14X,' S E',9X,' LAMBDA',6X,
     1 ' ESTIMATED PARTIALS USED ',/,5X,2E18.8,E13.3)
  380 FORMAT(/,' INCREMENTS ',5E18.8,/,(12X,5E18.8))
  381 FORMAT (13X,' PHI',10X,' LAMBDA',6X,' GAMMA ',6X,' LENGTH',/,
     1 5X,E18.8,3E13.3)
  382 FORMAT(1X,1E9.2,86X,1E9.2,/,1X,'+',99X,'+')
  383 FORMAT('1N = ',I3,5X,' K = ',I3,5X,'IP = ',I3,5X,' M = ',I3,5X,
     1 ' IFP = ',I3,5X,'GAMMA CRIT = ',E10.3,5X,'DEL = ',E10.3,/,
     2 ' FF = ',E10.3,5X,' T = ',E10.3,5X,' E = ',E10.3,5X,' TAU = ',
     3 E10.3,5X,' XL = ',E10.3,4X,'ZETA = ',E10.3,/)
  384 FORMAT(/,' PARAMETERS ',5E18.8,/,(12X,5E18.8))
  385 FORMAT (/,13X,' PHI',14X,' S E',9X,' LAMBDA',6X,
     1 ' ANALYTIC PARTIALS USED  ',/,5X,2E18.8,E13.3)
  386 FORMAT(//T12,'X(I)',T31,'Y OBS.',T49,'Y PRED.',T68,'DIFF',/)
  387 FORMAT (/,13X,' PHI',14X,' S E',11X,' LENGTH',6X,' GAMMA',6X,
     1 ' LAMBDA',6X,'ESTIMATED PARTIALS USED ',/,5X,2E18.8,3E13.3)
  388 FORMAT (/,13X,' PHI',14X,' S E',11X,' LENGTH',6X,' GAMMA',6X,
     1 ' LAMBDA',6X,'ANALYTIC PARTIALS USED  ',/,5X,2E18.8,3E13.3)
  389 FORMAT(2X,I3,' PARAMETER NOT USED ')
  390 FORMAT(2X,I3,' NONE FOUND ')
  391 FORMAT(2X,I3,36X,2E18.8  )
  392 FORMAT('  PTP INVERSE ')
  393 FORMAT(' ',/,' PARAMETER CORRELATION MATRIX ')
  394 FORMAT( 2X,I3,5E18.8)
  395 FORMAT(' ',/,' ',/,13X,' STD',17X,'ONE - PARAMETER ',21X,
     1 ' SUPPORT PLANE',/,3X,' B ',7X,' ERROR ',12X,' LOWER',12X,
     2 ' UPPER',12X,' LOWER',12X,' UPPER')
  396 FORMAT(' ',/,' ','  NONLINEAR CONFIDENCE LIMITS',/,/,
     1 ' PHI CRITICAL = ',E15.8)
  397 FORMAT(' ',/,'  PARA',6X,' LOWER B',8X,' LOWER PHI',10X,' UPPER B'
     1 ,8X,' UPPER PHI')
  398 FORMAT(' GAMMA LAMBDA TEST',5X,2E13.3)
  399 FORMAT(' EPSILON TEST ')
  400 FORMAT(' FORCE OFF ')
  401 FORMAT (5X,6E18.8,/,59X,2E18.8)
  402 FORMAT (' BAD DATA, SUBSCRIPTS FOR UNUSED BS = 0 ',/,/,/)
  403 FORMAT(2X,I3,5E18.8)
  404 FORMAT(' ',110A1)
  405 FORMAT(10A1)
  406 FORMAT (7F10.0)
  407 FORMAT (8F10.0)
  408 FORMAT('1')
  409 FORMAT('1N = ',I3,5X,' K = ',I3,5X,'IP = ',I3,5X,' M = ',I3,5X,
     1        /,' FF = ',E10.3,5X,' T = ',E10.3,
     2        5X,' E = ',E10.3,5X,' TAU = ',E10.3,/)
  410 FORMAT (' GAMMA EPSILON TEST   ')
  411 FORMAT (3X,I5,2X,10F10.4)
  412 FORMAT ('0 NEGATIVE DIAGONAL ELEMENT ')
      END
      BLOCK DATA
      COMMON/PLOT/IFP,PLOT(2)
      COMMON /PARAM/ PARAM(8)
      COMMON/WEITS/WTS(500000)
      DATA IFP,PLOT/0,2*0.0/
      DATA PARAM/4.,2.,.00005,.001,.01,45.,.00001,.1E-37/
      DATA WTS/500000*1.0/
      END
      SUBROUTINE CONLIM(N,K,M,Y,X,B,IP,IB,IBKS,IBD)
C     ..................................................................
C                  THIS IS THE CONFIDENCE LIMIT CALCULATION
      DIMENSION BS(50),SA(50),A(50,51)
      DIMENSION X(*),Y(*),B(*),IB(*)
       COMMON /PARAM/FF,T,E,TAU,AL,GAMCR,DEL,ZETA
       COMMON /AFCLSQ/ SE,PHI,PHIZ,WS,XL,IFSS2,IFSS3,IWS6,I,J,JJ,IBK2,
     ,IBKA,IBKM,IPLOT,BS,SA,A
C     TO INITIATE CONLIM
      IF (IBKS .EQ. 4)GO TO 21
      GO TO(15,17,18),IBKS
   15 IBKA1=IBKA-1
      GO TO(27,32),IBKA1
   17 GO TO 43
   18 IBK21=IBK2-1
      J=INDEX
      GO TO(158,27,125,134,144),IBK21
   21 DO 22 J=1,K
   22 B(J)=BS(J)
C     WRITE (26,201)N,K,IP,M,FF,T,E,TAU
      IBKA=2
C                  THIS WILL PRINT THE Y,YHAT,DELTA Y
      IBD=1
      GO TO 204
   27 IF (IPLOT .LE. 0) GO TO 32
   28 IBKA=3
      IPLOT=0
      IBD=1
      GO TO 204
   32 WS=N-K+IP
      SE=SQRT(PHI/WS)
      PHIZ=PHI
      IF (IFSS2.EQ.0)GO TO 38
C  36 WRITE (26,189)PHIZ,SE,XL
   36 CONTINUE
      GO TO 39
C  38 WRITE(26,190) PHIZ,SE,XL
   38 CONTINUE
C                  NOW WE HAVE MATRIX A
   39 CONTINUE
      IBKM=2
      IBD=2
      GO TO 204
C
C                  NOW WE HAVE C = A INVERSE
   43 DO 45 J=1,K
      IF(A(J,J).LT..0)GO TO 47
   45 SA(J)=SQRT(A(J,J))
      IBOUT=0
      GO TO 48
   47 IBOUT=1
   48 KST=-4
C     WRITE (26,194)
   50 KST=KST+5
      KEND=KST+4
      IF (KEND.LT.K)GO TO 54
      KEND=K
   54 DO 55 I=1,K
C  55 WRITE (26,196)I,(A(I,J),J=KST,KEND)
   55 CONTINUE
      IF (KEND.LT.K)GO TO 50
      IF (IBOUT.EQ.0)GO TO 61
C     WRITE (26,203)
      IBD=3
      GO TO 204
   61 DO 68 I=1,K
      DO 68 J=1,K
      WS=SA(I)*SA(J)
      IF(WS.GT. 0.)GOTO 67
   65 A(I,J)=0.
      GO TO 68
   67 A(I,J)=A(I,J)/WS
   68 CONTINUE
      DO 70 J=1,K
   70 A(J,J)=1.
C     WRITE (26,195)
      KST=-9
   73 KST=KST+10
      KEND=KST+9
      IF (KEND.LT.K)GO TO 77
      KEND=K
   77 DO 78 I=1,K
C  78 WRITE (26,202)I,(A(I,J),J=KST,KEND)
   78 CONTINUE
      IF (KEND.LT.K)GO TO 73
C                  GET T*SE*SQRT(C(I,I))
      DO 81 J=1,K
   81 SA(J)=  SE*SA(J)
C  82 WRITE (26,197)
   82 CONTINUE
      WS=K-IP
      DO 98 J=1,K
      IF (IP.LE.0)GO TO 89
   86 DO 88 I=1,IP
      IF (J.EQ.IB(I))GO TO 97
   88 CONTINUE
   89 HJTD=SQRT(WS*FF)*SA(J)
      STE=SA(J)
      OPL=BS(J)-SA(J)*T
      OPU=BS(J)+SA(J)*T
      SPL=BS(J)-HJTD
      SPU=BS(J)+HJTD
C     WRITE ( 26,200)J,STE,OPL,OPU,SPL,SPU
      GO TO 98
C  97 WRITE (26,191)J
   97 CONTINUE
   98 CONTINUE
C                  NONLINEAR CONFIDENCE LIMIT
      IF (IWS6.EQ.1) IBD=3
      IF (IWS6.EQ.1)GO TO 204
      WS=K-IP
      WS1=N-K+IP
      PKN=WS/WS1
      PC=PHIZ*(1.+FF*PKN)
C     WRITE (26,198)PC
C     WRITE (26,199)
      IFSS3=1
      J=1
  109 IBKP=1
      DO 112 JJ=1,K
  112 B(JJ)=BS(JJ)
      IF (IP.LE.0)GO TO 117
  114 DO 116 JJ=1,IP
      IF (J.EQ.IB(JJ))GO TO 173
  116 CONTINUE
  117 DD=-1.
      IBKN=1
  119 D=DD
      B(J)=BS(J)+D*SA(J)
      IBK2=4
      IBD=4
      INDEX=J
      GO TO 204
  125 PHI1=PHI
      IF (PHI1.GE.PC)GO TO 137
  127 D=D+DD
      IF (D/DD.GE.5.)GO TO 177
  129 B(J)=BS(J)+D*SA(J)
      IBK2=5
      IBD=4
      INDEX=J
      GO TO 204
  134 PHID=PHI
      IF (PHID.LT.PC)GO TO 127
      GO TO 146
  137 D=D/2.
      IF (D/DD.LE..001)GO TO 177
  139 B(J)=BS(J)+D*SA(J)
      IBK2=6
      IBD=4
      INDEX=J
      GO TO 204
  144 PHID=PHI
      IF (PHID.GT.PC)GO TO 137
  146 XK1=PHIZ/D+PHI1/(1.-D)+PHID/(D*(D-1.))
      XK2=-(PHIZ*(1.+D)/D+D/(1.-D)*PHI1+PHID/(D*(D-1.)))
      XK3=PHIZ-PC
      BC = (SQRT(XK2*XK2-4.*XK1*XK3)-XK2)/(2.*XK1)
      GO TO(151,153),IBKN
  151 B(J)=BS(J)-SA(J)*BC
      GO TO 154
  153 B(J)=BS(J)+SA(J)*BC
  154 IBK2=2
      IBD=4
      INDEX=J
      GO TO 204
  158 GO TO(159,164),IBKN
  159 IBKN=2
      DD=1.
      BL=B(J)
      PL=PHI
      GO TO 119
  164 BU=B(J)
      PU=PHI
      GO TO(167,169,171,175),IBKP
C 167 WRITE (26,196) J, BL, PL, BU, PU
  167 CONTINUE
      GO TO 185
C 169 WRITE (26,193) J, BU, PU
  169 CONTINUE
      GO TO 185
C 171 WRITE (26,196)J,BL, PL
  171 CONTINUE
      GO TO 185
C 173 WRITE (26,191)J
  173 CONTINUE
      GO TO 185
C 175 WRITE (26,192)J
  175 CONTINUE
      GO TO 185
  177 GO TO(178,180),IBKN
C                  DELETE LOWER PRINT
  178 IBKP=2
      GO TO 158
  180 GO TO(181,183),IBKP
C                  DELETE UPPER PRINT
  181 IBKP=3
      GO TO 158
C                  LOWER IS ALREADY DELETED, SO DELETE BOTH
  183 IBKP=4
      GO TO 158
  185 J=J+1
      J1=J-1
      IF (J1 .NE. K)GO TO 109
      DO 184 JJ=1,K
  184 B(JJ)=BS(JJ)
      IBD=5
  189 FORMAT (/,13X,' PHI',14X,' S E',9X,' LAMBDA',6X,
     1 ' ESTIMATED PARTIALS USED ',/,5X,2E18.8,E13.3)
  190 FORMAT (/,13X,' PHI',14X,' S E',9X,' LAMBDA',6X,
     1 ' ANALYTIC PARTIALS USED   ',/,5X,2E18.8,E13.3)
  191 FORMAT(2X,I3,' PARAMETER NOT USED ')
  192 FORMAT(2X,I3,' NONE FOUND ')
  193 FORMAT(2X,I3,36X,2E18.8  )
  194 FORMAT(' ',/,' PTP INVERSE ')
  195 FORMAT(' ',/,' PARAMETER CORRELATION MATRIX ')
  196 FORMAT(2X,I3,5E18.8)
  197 FORMAT(' ',/,' ',/,13X,' STD ',17X,' ONE - PARAMETER ',21X,
     1 ' SUPPORT PLANE',/,3X,' B',7X,' ERROR',12X,' LOWER',12X,
     2 ' UPPER',12X,' LOWER',12X,' UPPER')
  198 FORMAT(' ',/' ',/,'  NONLINEAR CONFIDENCE LIMITS ',/,/,
     1 ' PHI CRITICAL = ',E15.8)
  199 FORMAT(' ',/,'  PARA',6X,' LOWER B',8X,' LOWER PHI',10X,' UPPER B'
     1 ,8X,' UPPER PHI')
  200 FORMAT(2X,I3,5E18.8)
  201 FORMAT('1N = ',I3,5X,' K = ',I3,5X,'IP = ',I3,5X,' M = ',I3,5X,
     1       /,' FF = ',E10.3,5X,' T = ',E10.3,
     2       5X,' E = ',E10.3,5X,' TAU = ',E10.3,/)
  202 FORMAT (3X,I5,2X,10F10.4)
  203 FORMAT ('0 NEGATIVE DIAGONAL ELEMENT ')
  204 RETURN
      END
      SUBROUTINE GJR(A,N,EPS,MSING)
C GJR   OBJECT DECK DATE 9-18-73
C     GAUSS-JORDAN-RUTISHAUSER MATRIX INVERSION WITH DOUBLE PIVOTING.
      DIMENSION A(50,50),B(50),C(50),P(50),Q(50)
      INTEGER P,Q
      MSING=1
      DO 39 K=1,N
C     DETERMINATION OF THE PIVOT ELEMENT
      PIVOT=0.
      DO 13 I=K,N
      DO 13 J=K,N
      IF(ABS(A(I,J))-ABS(PIVOT))13,13,10
   10 PIVOT=A(I,J)
      P(K)=I
      Q(K)=J
   13 CONTINUE
      IF(ABS(PIVOT)-EPS)56,56,15
C     EXCHANGE OF THE PIVOTAL ROW WITH THE KTH ROW
   15 IF(P(K)-K)16,21,16
   16 DO 20 J=1,N
      L=P(K)
      Z=A(L,J)
      A(L,J)=A(K,J)
   20 A(K,J)=Z
C     EXCHANGE OF THE PIVOTAL COLUMN WITH THE KTH COLUMN
   21 IF(Q(K)-K)22,27,22
   22 DO 26 I=1,N
      L=Q(K)
      Z=A(I,L)
      A(I,L)=A(I,K)
   26 A(I,K)=Z
   27 CONTINUE
C     JORDAN STEP
      DO 36 J=1,N
      IF(J-K)33,30,33
   30 B(J)=1./PIVOT
      C(J)=1.
      GO TO 35
   33 B(J)=-A(K,J)/PIVOT
      C(J)=A(J,K)
   35 A(K,J)=0.
   36 A(J,K)=0.
      DO 39 I=1,N
      DO 39 J=1,N
   39 A(I,J)=A(I,J)+C(I)*B(J)
C     REORDERING THE MATRIX
      DO 54 M=1,N
      K=N-M+1
      IF(P(K)-K)43,48,43
   43 DO 47 I=1,N
      L=P(K)
      Z=A(I,L)
      A(I,L)=A(I,K)
   47 A(I,K)=Z
   48 IF(Q(K)-K)49,54,49
   49 DO 53 J=1,N
      L=Q(K)
      Z=A(L,J)
      A(L,J)=A(K,J)
   53 A(K,J)=Z
   54 CONTINUE
   55 RETURN
C  56 WRITE (26,57) P(K),Q(K),PIVOT
   56 CONTINUE
   57 format(' 0singular matrix   i=',i3,'j=',i3,'pivot=',e16.8/)
      MSING=2
      GO TO 55
      END

       subroutine fcode(y,x,b,f,i,N)
      dimension y(*), x(*), b(*)
      integer j
      real*4 vm(6)      
      call sdr2vm(b,vm)
      f = 0.0
      do j=1,6
         f = f + x(i+(j-1)*N)*vm(j)
      end do      
      return
      end

      subroutine sdr2vm(sdrM0,vm)
      dimension sdrM0(4),vm(6)
      real*4 ssi,sco,ssi2,sco2,dsi,dco,dsi2,dco2,rsi,rco
      real*4 DEG2RAD      
      DEG2RAD = 1.D0*datan(1.D0)/45.0
      ssi  = sin(sdrM0(1)*DEG2RAD)  
      sco  = cos(sdrM0(1)*DEG2RAD)  
      ssi2 = sin(2*sdrM0(1)*DEG2RAD)
      sco2 = cos(2*sdrM0(1)*DEG2RAD)
      dsi  = sin(sdrM0(2)*DEG2RAD)  
      dco  = cos(sdrM0(2)*DEG2RAD)  
      dsi2 = sin(2*sdrM0(2)*DEG2RAD)
      dco2 = cos(2*sdrM0(2)*DEG2RAD)
      rsi  = sin(sdrM0(3)*DEG2RAD)  
      rco  = cos(sdrM0(3)*DEG2RAD)  
      vm(1) =  sdrM0(4)*dsi2*rsi                           
      vm(2) = -sdrM0(4)*(dsi*rco*ssi2+dsi2*rsi*ssi*ssi)    
      vm(3) =  sdrM0(4)*(dsi*rco*ssi2-dsi2*rsi*sco*sco)    
      vm(4) = -sdrM0(4)*(dco*rco*sco +dco2*rsi*ssi)        
      vm(5) =  sdrM0(4)*(dco*rco*ssi -dco2*rsi*sco)        
      vm(6) = -sdrM0(4)*(dsi*rco*sco2+0.5 * dsi2*rsi*ssi2)   
      return
      end

      subroutine pcode(p,x,b,f,i)
      dimension p(*), x(*), b(*)
      return
      end

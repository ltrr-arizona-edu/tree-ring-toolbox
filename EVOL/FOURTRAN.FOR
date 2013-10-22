*RTN DFTPER: Discreet Fourier Transform of a series
*
* FORTRAN-77 code by R L Holmes, 10 JUL 1995
* Modified by R L Holmes, 20 JUL 1995
*
* Computes the discreet Fourier transform of a series, and its periodogram.
* The program is controlled by the following variables which are read in:
*
*   PT  Total proportion of the data to be tapered
*   KP  Degree of the polynomial to be removed, 0 or 1.
*
      PROGRAM DFTPER
      PARAMETER(MX=4096)
      CHARACTER QHD*40
      REAL X(MX),Y(MX)
*
* User input
      CALL QFPINI(KP,PT)
*
      DO I=1,MX
      X(I)=0.
      Y(I)=0.
      ENDDO ! I loop
*
* Read time series
      CALL DATIN(X,N,START,1.,QHD,7)
      JYR=NINT(START)
      WRITE(6,'('' Time span'',3I6,'' years'')')JYR,JYR+N-1,N
*
      CALL DETRND(X,N,KP)
      CALL TAPER(X,N,PT)
*
* NP2 is power of 2 .GE. N
      NP2=1
   20 NP2=NP2*2
      IF(NP2.LT.N)GOTO 20
*
      INV=0
      CALL FT01A(X,Y,NP2,INV)
*
      NPGM=NP2/2+1
      WRITE(QHD,'(''Fourier Real'')')
      CALL DATOUT(X,NPGM,START,1.,QHD,8)
      WRITE(QHD,'(''Fourier Imag'')')
      CALL DATOUT(Y,NPGM,START,1.,QHD,8)
      CON=(2.*FLOAT(NP2)/FLOAT(N))**2
*
      DO I=1,NPGM
      X(I)=(X(I)**2+Y(I)**2)*CON
      ENDDO ! I loop
*
      WRITE(QHD,'(''R-tilde'')')
      CALL DATOUT(X,NPGM,START,1.,QHD,9)
      STOP'DFTPER'
      END ! DFTPER
*RTN QFPINI: Keyboard input, opens files
*
* Programmed copyright (C) 06 JUL 1995 by
* Richard L Holmes, University of Arizona, Tucson, Arizona, USA
* Modified copyright (C) 20 JUL 1995
*
      SUBROUTINE QFPINI(KP,PT)
      CHARACTER CHAR,QF*12,QID*6
      WRITE(*,'(1X,A,''[0;1;32;40m'',A,''[2J'',A,
     &''[0;4H[1]  PROGRAM DFTPER'',T65,''Version TEST''/
     &'' Discreet Fourier Transform''/)')(CHAR(27),I=1,3)
      WRITE(*,'('' Run identification (to 6 letters or numbers): '',$)')
      READ(*,'(A)')QID
      LD=LARGO(QID)
*
* Open old file 7
      WRITE(*,'('' File with data series to be transformed: '',$)')
      READ(*,'(A)')QF
      OPEN(7,FILE=QF,STATUS='OLD',ERR=8)
      OPEN(6,FILE=QID(:LD)//'FT.OUT',STATUS='UNKNOWN')
      OPEN(8,FILE=QID(:LD)//'FT.DAT',STATUS='UNKNOWN')
      OPEN(9,FILE=QID(:LD)//'FT.RTL',STATUS='UNKNOWN')
*
    2 WRITE(*,'('' Proportion of series to be tapered  (0 to 1): '',$)')
      READ(*,*,ERR=2)PT
      IF(PT.LT.0..OR.PT.GT.1.)GOTO 2
    3 WRITE(*,'('' Degree of polynomial to be removed  (0 or 1): '',$)')
      READ(*,*,ERR=3)KP
      KP=MAX0(0,MIN0(1,KP))
*
      WRITE(6,'(
     &'' PROGRAM DFTPER - Discreet Fourier Transform and Periodogram''/
     &1X,A,'' Data file''/'' Proportion of data tapered'',F8.4/
     &'' Least-squares polynomial of degree'',I2,'' removed'')')QF,PT,KP
      IF(KP.EQ.0)THEN
      WRITE(6,'('' (Mean subtracted)'')')
      ELSE
      WRITE(6,'('' (Least-squares regression line subtracted)'')')
      ENDIF
*
      WRITE(*,'(/'' Output files are '',A,''FT.OUT''/
     &18X,A,''FT.DAT''/14X,''and '',A,''FT.RTL'')')(QID(:LD),K=1,3)
      WRITE(6,'(/'' Output files are '',A,''FT.OUT, '',A,
     &''FT.DAT and '',A,''FT.RTL'')')(QID(:LD),K=1,3)
      RETURN
    8 WRITE(*,'(/'' File: '',A)')QF
      STOP'  cannot be opened'
      END ! QFPINI
*RTN DATIN: Data input
*
* FORTRAN-77 code by R L Holmes, 06 JUL 1995
* Modified by R L Holmes, 10 JUL 1995
*
* Arguments:
*
*   X     Time series
*   N     Length of series
*   START Time value of first data point
*   STEP  Time increment between data points
*   IO    Logical unit for input data file
*
      SUBROUTINE DATIN(X,N,START,STEP,QHD,IO)
      PARAMETER(MX=4096)
      CHARACTER QHD*40,QFMT*72
      REAL X(MX)
      QFMT='Compact'
      CALL DATR(QHD,N,X,JYR,IO)
      IF(N.LE.0)RETURN
      START=FLOAT(JYR)
      END=START+FLOAT(N-1)*STEP
      LH=LARGO(QHD)
      LF=LARGO(QFMT)
*
      WRITE(6,'(/'' Header '',A/'' Series time span'',3I6,'' years'',
     &''  Data format '',A/'' Data time span'',2F10.3,''  Increment'',
     &F7.3)')QHD(:LH),JYR,JYR+N-1,N,QFMT(:LF),START,END,STEP
      END ! DATIN
*RTN DATOUT: Data output
*
* FORTRAN-77 code by R L Holmes, 06 JUL 1995
* Modified by R L Holmes, 10 JUL 1995
*
* Outputs series X in the format expected by Routine DATIN.
* The Header line is blank except for sequencing.
*
* Arguments
*
*   X      Time series
*   N      Length of series
*   START  Origin of time series
*   STEP   Time series increment
*   IO     Logical unit for output
*
      SUBROUTINE DATOUT(X,N,START,STEP,QHD,IO)
      CHARACTER QHD*40
      REAL X(N)
*##      WRITE(IO,'(T80,''1''/I5,T80,''2''/'' (5E15.8)'',T80,''3''/
*##     &2F10.5,T80,''4'')')N,START,STEP
*##      K=4
*##      L=5
*##      IHI=0
*##   10 ILO=IHI+1
*##      IHI=IHI+L
*##      K=K+1
*##      IF(IHI.GT.N)GOTO 20
*##      WRITE(IO,'(5E15.8,I5)')(X(I),I=ILO,IHI),K
*##      GOTO 10
*##   20 IF(ILO.GT.N)RETURN
*##      LIM=IHI-N
*##      WRITE(IO,'(5E15.8,I5)')(X(I),I=ILO,N),(Z(I),I=1,LIM),K
*
      WRITE(QHD(13:),'(F10.3,'' Start'',F7.3,'' Step'')')START,STEP
      JYR=NINT(START)
      CALL DATW(QHD,N,X,JYR,IO)
      END ! DATOUT
*RTN DETRND: Subtracts mean or linear regression from series
*
* FORTRAN-77 code by R L Holmes, 11 JUL 1995
* Modified by R L Holmes, 11 JUL 1995
*
* Computes and subtracts from the series X either the series mean or the
* least-squares straight line.  Note that these are the least-squares
* polynomials of degree 0 and 1 respectively.
*
* Arguments:
*  X   Time series
*  N   Length of series
*  KP  Degree of polynomial to be subtracted from the series:
*        0 = Mean subtracted
*        1 = Least-squares regression line subtracted
*
      SUBROUTINE DETRND(X,N,KP)
      REAL X(N)
*
* Remove mean
      SUMX=0.
      DO I=1,N
      SUMX=SUMX+X(I)
      ENDDO ! I loop
      XBAR=SUMX/FLOAT(N)
*
      DO I=1,N
      X(I)=X(I)-XBAR
      ENDDO ! I loop
      IF(KP.LE.0)RETURN
*
* Remove trend
      TBAR=FLOAT(N+1)/2.
      SUMTT=FLOAT(N*(N*N-1))/12.
*
      SUMTX=0.
      DO I=1,N
      SUMTX=SUMTX+X(I)*(FLOAT(I)-TBAR)
      ENDDO ! I loop
      BETA=SUMTX/SUMTT
*
      DO I=1,N
      X(I)=X(I)-BETA*(FLOAT(I)-TBAR)
      ENDDO ! I loop
      END ! DETRND
*RTN TAPER: Applies split-cosine-bell tapering to the time series
*
* FORTRAN-77 code by R L Holmes, 11 JUL 1995
* Modified by R L Holmes, 11 JUL 1995
*
* Applies split-cosine-bell tapering to the time series X.
* The argument P is the total proportion of the series which is tapered.
*
* Arguments
*
*    X   Time series, tapered on exit
*    N   Length of series
*    PT  Total proportion of the series to be tapered
*
      SUBROUTINE TAPER(X,N,PT)
      REAL X(N)
      DATA PI/3.14159265/
      IF(PT.LE.0..OR.PT.GT.1.)RETURN
      M=NINT(PT*FLOAT(N))/2
      DO I=1,M
      WT=.5-.5*COS(PI*(FLOAT(I)-.5)/FLOAT(M))
      X(I)=X(I)*WT
      X(N+1-I)=X(N+1-I)*WT
      ENDDO ! I loop
      END ! TAPER
*RTN FT01A: Implements the Sande-Tukey radix-2 fast Fourier transform
*
* FORTRAN-77 code by R L Holmes, 11 JUL 1995
* Modified by R L Holmes, 21 JUL 1995
*
* Implements the Sande-Tukey radix-2 fast Fourier transform.
* Either the direct or the inverse transform may be computed.
*
* Arguments
*
*  XR  Real part of the series and the transform
*  XI  Imaginary part of the series and the transform
*  N   Length of series
*  INV  = 0 for direct transform
*       = 1 for inverse transform
*
* The direct transform is:
*                  N
*         (1/N) * SUM X(T)*EXP(-2*PI*I*(T-1)*(J-1)/N), J=1,N
*                 T=1
*
* The inverse transform is:
*                  N
*                 SUM X(T)*EXP( 2*PI*I*(T-1)*(J-1)/N), J=1,N
*                 T=1
*
      SUBROUTINE FT01A(XR,XI,N,INV)
      REAL XR(N),XI(N),UR(15),UI(15)
*
* Define UR and UI arrays
      UR(1)=0.
      UI(1)=1.
*
      DO I=2,15
      UR(I)=SQRT((1.+UR(I-1))/2.)
      UI(I)=UI(I-1)/(2.*UR(I))
      ENDDO ! I loop
*
      IF(N.LE.0.OR.N.GT.2**16)THEN
      WRITE(6,'('' >> RTN FT01A: N out of range'',I12)')N
      STOP'>> Error'
      ENDIF
*
      N0=1
      II=0
  140 N0=N0+N0
      II=II+1
      IF(N0.LT.N)GOTO 140
*
      I1=N0/2
      I3=1
      I0=II
*
      DO 260 I4=1,II
*
      DO 250 K=1,I1
      WR=1.
      WI=0.
      KK=K-1
*
      DO 230 I=1,I0
      IF(KK.EQ.0)GOTO 240
      IF(MOD(KK,2).EQ.0)GOTO 230
      J0=I0-I
      WS=WR*UR(J0)-WI*UI(J0)
      WI=WR*UI(J0)-WI*UR(J0)
      WR=WS
  230 KK=KK/2
*
  240 IF(INV.EQ.0)WI=-WI
      L=K
*
      DO 250 J=1,I3
      L1=L+I1
      ZR=XR(L)+XR(L1)
      ZI=XI(L)+XI(L1)
      Z     =WR*(XR(L)-XR(L1))-WI*(XI(L)-XI(L1))
      XI(L1)=WR*(XI(L)-XI(L1))+WI*(XR(L)-XR(L1))
      XR(L1)=Z
      XR(L)=ZR
      XI(L)=ZI
  250 L=L1+I1
      I0=I0-1
      I3=I3+I3
  260 I1=I1/2
      UM=1.
      IF(INV.EQ.0)UM=1./FLOAT(N0)
*
      DO 310 J=1,N0
      K=0
      J1=J-1
*
      DO 320 I=1,II
      K=2*K+MOD(J1,2)
  320 J1=J1/2
*
      K=K+1
      IF(K.LT.J)GOTO 310
      ZR=XR(J)
      ZI=XI(J)
      XR(J)=XR(K)*UM
      XI(J)=XI(K)*UM
      XR(K)=ZR*UM
      XI(K)=ZI*UM
  310 CONTINUE
*
      END ! FT01A
*RTN DATR: Read time series in compact format (~ or @)
*
* Programmed copyright (C) FEB 1990 by
* Richard L Holmes, University of Arizona, Tucson, Arizona, USA
* Modified copyright (C) 30 AUG 1994
*
* ID:  Identification up to 40 letters
* N:   Number of data in series; N returned as negative at end of file
* Y:   Real time series
* JYR: First year of time series
* IO:  Logical unit for reading data
*
      SUBROUTINE DATR(ID,N,Y,JYR,IO)
      PARAMETER(MX=4096)
      IMPLICIT CHARACTER(Q)
      COMMON/DPLQ/QFS(64),QFF(64),QQ
      COMMON/DPLN/IU(64),IUO,IAJ
      CHARACTER QQ*254,LIN*80,FRM*8,ID*(*)
      REAL Y(1)
      DATA QT/' '/
*
* Check format
      N=0
      LZ=MIN0(LEN(ID),40)
    2 READ(IO,'(A)',ERR=99,END=99)LIN
      LL=LARGO(LIN)
*
      IF(LL.LE.1)THEN
      GOTO 2
      ELSEIF(LIN(LL:LL).EQ.'~')THEN
      GOTO 1
      ELSEIF(LIN(LL:LL).EQ.'@')THEN
      FRM=LIN(LL-7:LL-5)//'F'//LIN(LL-2:LL-2)//'.0)'
      LIN(LL-11:LL-1)='  0'//FRM
      GOTO 1
      ELSE
      IF(QQ(14:14).EQ.' ')WRITE(*,'('' Title:''/1X,A)')LIN
      IF(IUO.GT.0.AND.QQ(14:14).EQ.' ')THEN
      WRITE(IUO,'(/'' Title of data: '',A)')LIN
      QT='1'
      ENDIF
      GOTO 2
      ENDIF
*
* Read header and data
    1 READ(LIN,'(I8,T11,I8)',ERR=2,END=99)N,JYR
      READ(LIN(LL-11:LL-1),'(I3,A)',ERR=2,END=99)NC,FRM
      LIN(LL-11:LL)=' '
      ID(1:LZ)=LIN(22:LZ+21)
      READ(IO,FRM,ERR=97,END=98)(Y(J),J=1,N)
      FAC=10.**NC
*
      DO J=1,N
      Y(J)=Y(J)*FAC
      ENDDO
      QQ(14:14)=QT
      RETURN
*
* Error messages
   97 WRITE(*,'('' >> Sbr DATR:  Read error'')')
      RETURN
   98 WRITE(*,'('' >> Sbr DATR:  EOF unexpected'')')
      RETURN
   99 N=-1
      QQ(14:14)=QT
      END ! DATR
*RTN DATW: Write time series on disk file in compact format (~)
*
* Programmed copyright (C) FEB 1990 by
* Richard L Holmes, University of Arizona, Tucson, Arizona, USA
* Modified copyright (C) 11 FEB 1994
*
* ID:  Series ident to 40 letters
* N:   Number of data in series
* Y:   Real time series
* JYR: First year of time series
* IO:  Logical unit for writing data
*
      SUBROUTINE DATW(ID,N,Y,JYR,IO)
      PARAMETER(MX=4096)
      COMMON/IARR/IZZ(MX)
      CHARACTER FRM*8,ID*(*)
      REAL Y(MX)
      IF(N.LE.0)RETURN
      LA=MIN0(LEN(ID),40)
      LZ=LARGO(ID(1:LA))
*
* Encontrar espacio necesario:
      NE=1
      LMAX=1
      VMAX=0.
      FAC=1.
*
      DO I=1,N
      FAC=AMAX1(FAC,Y(I),Y(I)*(-10.))
      IF(ABS(Y(I)).GT.VMAX)THEN
      VMAX=ABS(Y(I))
      LMAX=I
      ENDIF
      ENDDO ! I loop
*
* Factor para multiplicar datos
      FAC1=1.E16
      DO K=1,32
      FAC1=FAC1*.1
      IF(FAC*FAC1.LT.999999.5)GOTO 2
      ENDDO ! K loop
*
* Convertir arreglo a entero por factor
    2 DO J=1,N
      IZZ(J)=NINT(Y(J)*FAC1)
      ENDDO ! J loop
*
* Calcular acortamiento - incrementar modulo hasta hallar discrepancia
      M=1
      DO I=1,6
      M=M*10
      DO J=1,N
      IF(MOD(IZZ(J),M).NE.0)THEN
      M=M/10
      NE=7-I
      GOTO 1
      ENDIF
      ENDDO ! J loop
      ENDDO ! I loop
*
* Convertir datos enteros
    1 IF(M.GE.10)THEN
      DO J=1,N
      IZZ(J)=IZZ(J)/M
      ENDDO ! J loop
      ENDIF
*
      IF(IZZ(LMAX).NE.0)THEN
      V=ABS(VMAX/FLOAT(IZZ(LMAX)))
      ELSE
      V=0.
      ENDIF
*
      IF(V.GT.0.)THEN
      NC=NINT(ALOG10(V))
      ELSE
      NC=0
      ENDIF
*
* Primera linea:  N (1:8), JYR (11:18), ID(22:61),
* Exponente (69:71), Formato (72:79), Clave (~) (80)
      NL=80/NE
      WRITE(FRM,'(''('',I2,''F'',I1,''.0)'')')NL,NE
      WRITE(IO,'(I8,''=N'',I8,''=I'',T22,A,T69,I3,A,''~'')')
     &N,JYR,ID(1:LZ),NC,FRM
      WRITE(FRM,'(''('',I2,''I'',I1,'')  '')')NL,NE
*
* Escribir serie de datos
      WRITE(IO,FRM)(IZZ(J),J=1,N)
      END ! DATW
*RTN LARGO: Position of last non-blank character in string
*
* Programmed copyright (C) FEB 1991 by
* Richard L Holmes, University of Arizona, Tucson, Arizona, USA
* Modified copyright (C) 01 DEC 1992
*
      FUNCTION LARGO(Q)
      CHARACTER Q*(*)
      LL=LEN(Q)
      DO J=LL,1,-1
      K=ICHAR(Q(J:J))
      IF(K.NE.32.AND.K.NE.0)THEN
      LARGO=J
      RETURN
      ENDIF
      ENDDO ! J loop
      LARGO=1
      END ! LARGO

C   ANTI-ROBINSON SERIATION
C   branch-and-bound
C   by Brusco, and Stahl, S.
C   R Interface by Michael Hahsler


C      PROGRAM DYNAMIC
C      SUBROUTINE dynamic(N, A, EPS, X)
      SUBROUTINE bbwrcg(N, A, EPS, X, Q, D, DD, S, UNSEL, IVERB)
      IMPLICIT INTEGER(A-Z)
C      DOUBLE PRECISION TIMEA,TIMEB,TIMTOT,A(50,50),EPS
      DOUBLE PRECISION EPS, A(N,N), D(N,N,N),
     1    DD(N,N,N),ZBEST,Z,ACT,DELTA,ZBD,IDX1,IDX2
      REAL S1
      INTEGER X(N),Q(N),S(N),UNSEL(N)

C     EPS is unused this is to supress the warning.
      EPS = 1.0d-07

C       Initialize R RNG
      CALL getrngstate()

      OLDM=0
      CHECKS=0


C
C #################################################################
C 10/13/01 This program fits an "weighted" row gradient criterion
C         to a symmetric proximity matrix.  Count +1 if the anti-
C         Robinson triple is satisfied, -1 if its not, and 0 for
C         ties.  Only look at upper half of matrix
C 07/20/02: Improved symmetry test implemented.
C 07/26/03: Fixed the incorrect symmetry test, added an interchange test
C          avoid use of so many "IF" statements using F & D matrices
C 12/24/03: Add insertion test to interchange test.
C 07/09/15: Fixed memory issue (MFH)
C #################################################################
C
C      OPEN(1,FILE='AMAT.DAT')          ! Dissimilarity matrix
C      OPEN(2,FILE='SEQ.OUT')           ! Output file
C      EPS = 1.0d-07
C      READ(1,*) N                      ! Read number of objects
C      WRITE(*,*) 'TYPE 1 FOR HALF MATRIX OR TYPE 2 FOR FULL MATRIX'
C      READ(*,*) ITYPE
C      ITYPE = 2
C      IF(ITYPE.EQ.2) THEN
C        READ(1,*) ((A(I,J),J=1,N),I=1,N)
C      ELSE
C        DO J = 2,N
C          READ(1,*) (A(I,J),I=1,J-1)
C        END DO
C        DO J = 2,N
C          DO I = 1,J-1
C            A(J,I) = A(I,J)
C          END DO
C        END DO
C      END IF
C      CALL GETTIM (IHR, IMIN, ISEC, I100)
C      CALL GETDAT (IYR, IMON, IDAY)
C      TIMEA=DFLOAT(86400*IDAY+3600*IHR+60*IMIN+ISEC)+DFLOAT(I100)/100.
      DO I = 1,N
        A(I,I) = 0.0D0
      END DO
C
      DO 848 I = 1,N
        DO 849 J = 1,N
          IF(I.EQ.J) GO TO 849
          DO 850 K = 1,N
            IF(I.EQ.K.OR.J.EQ.K) GO TO 850
C bbwrg
C            D(I,J,K) = A(I,K) - A(I,J)
            D(I,J,K) = 2.*A(I,K) - A(I,J) - A(J,K)


 850      CONTINUE
 849    CONTINUE
 848  CONTINUE
C
      DO 851 I = 1,N
        DO 852 J = 1,N
          IF(I.EQ.J) GO TO 852
          DO 853 K = 1,N
            IF(I.EQ.K.OR.J.EQ.K) GO TO 853
            ACT=D(I,J,K)
            IF(D(I,K,J).GT.ACT) ACT = D(I,K,J)
            IF(D(J,I,K).GT.ACT) ACT = D(J,I,K)
            DD(I,J,K) = ACT
 853      CONTINUE
 852    CONTINUE
 851  CONTINUE
C Run heuristic to find a good objective value
      IF (IVERB == 1) THEN
        CALL intpr('Run heuristic', -1, IVERB, 0)
      ENDIF

      ZBEST = 0.0D0
C      DO 3500 JJJ = 1,100
      DO 3500 JJJ = 1,100
        DO I = 1,N
          UNSEL(I) = I
          Q(I) = 0
        END DO
        NNSEL = N
C 3501   CALL RANDOM(S1)
C 3501   S1 = rand()
 3501   CALL unifrand(S1)
        ISEL = INT(1. + S1*FLOAT(NNSEL))
        IF(ISEL.GT.NNSEL) ISEL = NNSEL
        Q(NNSEL) = UNSEL(ISEL)
        DO J = ISEL,NNSEL-1
          UNSEL(J) = UNSEL(J+1)
        END DO
        NNSEL = NNSEL - 1
        IF(NNSEL.GT.0) GO TO 3501
C        WRITE(*,72) (Q(J),J=1,N)
C 72     FORMAT(20I3)
        Z = 0.0D0
        DO I = 1,N-2
          R1 = Q(I)
          DO J = I+1,N-1
            R2 = Q(J)
            DO K = J+1,N
              R3 = Q(K)
              Z = Z + D(R1,R2,R3)
            END DO
          END DO
        END DO
 3502   ITRIG = 0
        DO II = 1,N-1
          DO JJ = II+1,N
C   R interrupt
            CALL rchkusr()
C
            R3 = Q(JJ)
            R2 = Q(II)
            DELTA=0.0D0
            DO I = 1,II-1
              R1 = Q(I)
              DELTA = DELTA + D(R1,R3,R2) - D(R1,R2,R3)
              DO J = II+1,JJ-1
                R4 = Q(J)
                DELTA = DELTA + D(R1,R3,R4) - D(R1,R2,R4)
                DELTA = DELTA + D(R1,R4,R2) - D(R1,R4,R3)
              END DO
            END DO
            DO J = II+1,JJ-1
              R4 = Q(J)
              DELTA = DELTA + D(R3,R4,R2) - D(R2,R4,R3)
              DO K = JJ+1,N
                R5 = Q(K)
                DELTA = DELTA + D(R4,R2,R5) - D(R4,R3,R5)
                DELTA = DELTA + D(R3,R4,R5) - D(R2,R4,R5)
              END DO
            END DO
            DO K = JJ + 1,N
              R5 = Q(K)
              DELTA = DELTA + D(R3,R2,R5) - D(R2,R3,R5)
            END DO
            DO I = II+1,JJ-2
              DO J = I+1,JJ-1
                R4A = Q(I)
                R4B = Q(J)
                DELTA = DELTA + D(R4A,R4B,R2) - D(R4A,R4B,R3)
                DELTA = DELTA + D(R3,R4A,R4B) - D(R2,R4A,R4B)
              END DO
            END DO
            IF(DELTA.GT.0) THEN
              Z = Z + DELTA
              Q(II) = R3
              Q(JJ) = R2
              ITRIG = 1
            END IF
          END DO
        END DO
        IF(ITRIG.EQ.1) GO TO 3502
        IF(Z.GT.ZBEST) ZBEST = Z

 3500 CONTINUE
C      WRITE(2,3505) ZBEST
      IF (IVERB == 1) THEN
C          WRITE(*,3505) ZBEST
          CALL dblepr('HEURISTIC OBJ VALUE', -1, DBLE(ZBEST), 1)
      ENDIF
C 3505 FORMAT(' HEURISTIC OBJ VALUE ',F20.4)
      Z = ZBEST-1
      DO I = 1,N
        Q(I) = 0
      END DO
C
      M=1
      Q(M)=1
      S(1)=1
      trig=1
      DO K = 2,N
        Q(K)=0
      END DO
C
  1   M = M + 1
C
C
      CHECKS=CHECKS+1
      IF (IVERB == 1 .AND. M .GT. OLDM) THEN
C          WRITE (*,6000) M+1, CHECKS
          CALL intpr('reached position', -1, M+1, 1)
          CALL intpr('with following number of checks', -1, CHECKS, 1)
C 6000 FORMAT('reached position ', I5, ' with ', I9, ' checks')
          OLDM=M
      ENDIF

C
C   R interrupt
      CALL rchkusr()
C
C   main loop
C
  2   Q(M)=Q(M)+1
C
C MFH: Make sure to not get out of bounds with S(Q(M)) - 7/9/15
      IF(Q(M).GT.N) GO TO 222
      IF(S(Q(M)).EQ.1) GO TO 2               ! REDUNDANCY
  222 IF(M.EQ.1.AND.Q(M).GT.N) GO TO 9       ! TERMINATE
      IF(M.GT.1.AND.Q(M).GT.N) GO TO 7       ! GO TO RETRACTION
C only for bbwrcg
      IF(TRIG.EQ.0.AND.Q(M).EQ.2) GO TO 2    ! SYMMETRY FATHOM
C
      S(Q(M))=1
      IF(M.EQ.1) GO TO 1
      IF(M.EQ.N-1) THEN
        CALL EVALBBWRCG(ZBD,Q,N,D)
        IF(ZBD.GT.Z) THEN
          Z=ZBD
          IF (IVERB == 1) THEN
C              WRITE(*,*) 'Eval =',z
              CALL dblepr('Eval', -1, DBLE(z), 1)
          ENDIF
          DO I = 1,N
            X(I)=Q(I)
          END DO
        END IF
        Q(N)=0
        S(Q(M))=0
        GO TO 2
      ELSE
        DO 251 MM = M-1,1,-1    ! Insertion Test
          R3=Q(M)
          IDX1=0
          IDX2=0
          DO I = 1,MM-1
            R1=Q(I)
            DO J = MM,M-1
              R4=Q(J)
              IDX1=IDX1+D(R1,R4,R3)
              IDX2=IDX2+D(R1,R3,R4)
C
            END DO
C
          END DO
C
          DO 250 I = 1,N
            IF(S(I).EQ.1) GO TO 250
            R5=I
C
            DO J = MM,M-1
              R4=Q(J)
              IDX1=IDX1+D(R4,R3,R5)
              IDX2=IDX2+D(R3,R4,R5)
            END DO
C
  250     CONTINUE
C
          DO J = MM, M-2
            DO JJ = J+1, M-1
              R4A=Q(J)
              R4B=Q(JJ)
              IDX1=IDX1+D(R4A,R4B,R3)
              IDX2=IDX2+D(R3,R4A,R4B)
            END DO
          END DO
          IF(IDX1.LT.IDX2) THEN
            S(Q(M))=0
C            ism2 = ism2 + 1
            GO TO 2
          END IF
  251   CONTINUE
C        go to 253
C
        DO 151 MM = M-2,1,-1    ! Interchange Test
          R3=Q(M)
          R2=Q(MM)
          IDX1=0
          IDX2=0
          DO J = MM+1,M-1
            R4 = Q(J)
            IDX1=IDX1+D(R2,R4,R3)
            IDX2=IDX2+D(R3,R4,R2)
          END DO
          DO I = 1,MM-1
            R1=Q(I)
            IDX1=IDX1+D(R1,R2,R3)
            IDX2=IDX2+D(R1,R3,R2)
            DO J = MM+1,M-1
              R4=Q(J)
              IDX1=IDX1+D(R1,R2,R4)
              IDX2=IDX2+D(R1,R3,R4)
C
              IDX1=IDX1+D(R1,R4,R3)
              IDX2=IDX2+D(R1,R4,R2)
C
            END DO
C
          END DO
C
          DO 150 I = 1,N
            IF(S(I).EQ.1) GO TO 150
            R5=I
            IDX1=IDX1+D(R2,R3,R5)
            IDX2=IDX2+D(R3,R2,R5)
C
            DO J = MM+1,M-1
              R4=Q(J)
              IDX1=IDX1+D(R2,R4,R5)
              IDX2=IDX2+D(R3,R4,R5)
C
              IDX1=IDX1+D(R4,R3,R5)
              IDX2=IDX2+D(R4,R2,R5)
            END DO
C
  150     CONTINUE
C
          DO J = MM+1, M-2
            DO JJ = J+1, M-1
              R4A=Q(J)
              R4B=Q(JJ)
              IDX1=IDX1+D(R4A,R4B,R3)
              IDX2=IDX2+D(R4A,R4B,R2)
              IDX1=IDX1+D(R2,R4A,R4B)
              IDX2=IDX2+D(R3,R4A,R4B)
            END DO
          END DO
          IF(IDX1.LT.IDX2) THEN
C            ism = ism + 1
            S(Q(M))=0
            GO TO 2
          END IF
  151   CONTINUE
C
        CALL BOUND2BBWRCG(ZBD,N,Q,M,D,S,DD)
        IF(ZBD.LE.Z) THEN
          S(Q(M))=0
C          ism3 = ism3 + 1
          GO TO 2
        END IF
        IF(Q(M).EQ.1) TRIG=1
        GO TO 1
      END IF
C
   7  IF(Q(M).EQ.1) TRIG=0
C MFH: Make sure to not get out of bounds with S(Q(M)) - 6/9/15
      IF(Q(M).GT.N) GO TO 777
      S(Q(M))=0
  777 Q(M)=0
      M=M-1
      IF(Q(M).EQ.1) TRIG=0
      S(Q(M))=0
C
C      WRITE(*,*) 'X',(X(J),J=1,N)
C      WRITE(*,*) 'Q',(Q(J),J=1,N)
C

      GO TO 2
C   9  CALL GETTIM (IHR, IMIN, ISEC, I100)
C      CALL GETDAT (IYR, IMON, IDAY)
C      TIMEB=DFLOAT(86400*IDAY+3600*IHR+60*IMIN+ISEC)+DFLOAT(I100)/100.
C      TIMTOT=TIMEB-TIMEA
C      write(*,*) ism,ism2,ism3
C      WRITE(*,69) Z,TIMTOT
C  9   WRITE(*,69) Z
C      WRITE(2,69) Z,TIMTOT
C      WRITE(2,70) (X(I),I=1,N)
C 69   FORMAT(' MAXIMUM UNWEIGHTED ROW GRADIENT INDEX ',I20)
C 69   FORMAT(' MAXIMUM UNWEIGHTED ROW GRADIENT INDEX ',I7,' CPU TIME ',
C     1        F8.2)
C 70   FORMAT(30I3)
C
  9   IF (IVERB == 1) THEN
C          PRINT *, 'total number of checks: ', CHECKS
          CALL intpr('total number of checks', -1, CHECKS, 1)
      ENDIF

C    Return R RNG
      CALL Putrngstate()


      RETURN
      END
C
      SUBROUTINE BOUND2BBWRCG(ZBD,N,Q,M,D,S,DD)
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION D(N,N,N),ZBD,DD(N,N,N),Z1,Z2,Z3,ZA,ZB,
     1   ZCT,N4
C       ACT is now unused

      INTEGER Q(N),S(N)
      Z1=0
      DO I = 1,M-2
        R1=Q(I)
        DO J = I+1,M-1
          R2=Q(J)
          DO K = J+1,M
            R3=Q(K)
            Z1=Z1+D(R1,R2,R3)
          END DO
        END DO
      END DO
C
      Z2=0
      DO I = 1,M-1
        R1=Q(I)
        DO J = I+1,M
          R2=Q(J)
          DO 60 K = 1,N
            IF(S(K).EQ.1) GO TO 60
            R3=K
            Z2=Z2+D(R1,R2,R3)
 60       CONTINUE
        END DO
      END DO
C
      Z3=0
      DO 90 I = 1,N-1
        IF(S(I).EQ.1) GO TO 90
        R2=I
        DO 91 J = I+1,N
          IF(S(J).EQ.1) GO TO 91
          R3=J
          ZA=0
          ZB=0
          DO 92 K = 1,M
            R1=Q(K)
            ZA=ZA+D(R1,R2,R3)
            ZB=ZB+D(R1,R3,R2)
 92       CONTINUE
          ZCT=ZA
          IF(ZB.GT.ZCT) ZCT=ZB
          Z3=Z3+ZCT
 91     CONTINUE
 90   CONTINUE
C
      N4=0
      DO 93 I = 1,N-2
        IF(S(I).EQ.1) GO TO 93
        R1=I
        DO 94 J = I+1,N-1
          IF(S(J).EQ.1) GO TO 94
          DO 95 K = J+1,N
            IF(S(K).EQ.1) GO TO 95
C            ACT=D(I,J,K)
C            IF(D(I,K,J).GT.ACT) ACT=D(I,K,J)
C            IF(D(J,I,K).GT.ACT) ACT=D(J,I,K)
C            N4=N4+ACT
             N4 = N4 + DD(I,J,K)
  95      CONTINUE
  94    CONTINUE
  93  CONTINUE
C      N1=N*(N-1)*(N-2)/3    ! This bound is OK!  The N1 is total
C      N2=M*(M-1)*(M-2)/3    ! and N2 and N3 are truly computed terms.
C      N3=(N-M)*(M*(M-1))  ! So N1-N2-N3 assumes +1 for rest, which
      ZBD=Z1+Z2+Z3+n4   ! (N-M)*(N-M-1)*(N-M-2)/3      +n4
C      WRITE(*,98) N,M,N1,N2,N3,Z1,Z2,N1-N2-N3,ZBD
C 98   FORMAT(9I7)
      RETURN
      END
C
      SUBROUTINE EVALBBWRCG(ZBD,Q,N,D)
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION D(N,N,N),ZBD
      INTEGER Q(N)
      ZBD=0
      DO 85 I = 1,N
        DO J = 1,N-1
          IF(Q(J).EQ.I) GO TO 85
        END DO
        Q(N)=I
 85   CONTINUE
      DO I = 1,N-2
        R1=Q(I)
        DO J = I+1,N-1
          R2=Q(J)
          DO K = J+1,N
            R3=Q(K)
            ZBD=ZBD+D(R1,R2,R3)
          END DO
        END DO
      END DO
      RETURN
      END

C   ANTI-ROBINSON SERIATION
C   simulated annealing algorithm - provides an initial permutation
C   by Brusco, M., Koehn, H.F., and Stahl, S.
C   R Interface by Michael Hahsler

C      PROGRAM SANNEAL
      SUBROUTINE arsa(N, A, COOL, TMIN, NREPS, IPERM, R1, R2, D, U, 
     1 S, T, SB, IVERB)
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N,N)
      DIMENSION IPERM(N)
C      DOUBLE PRECISION A(400,400), SOLS(100), RMED(100),
      DOUBLE PRECISION R1(N*N/2), R2(N*N/2), D(N,N)
      REAL S1, RCRIT
      INTEGER U(N), S(N), UNSEL, T(100,N), SB(N), Q

C	Initialize R RNG
      CALL getrngstate()

      IF (IVERB == 1) THEN
C         PRINT *, 'Anti-Robinson seriation by simulated annealing'
C         PRINT *, 'based on arsa.f by Brusco, M., Koehn, H.F.,',
C     1 'and Stahl, S. (2007)'

          CALL FPRINTF('Anti-Robinson seriation by simulated '
     1 //'annealing', 46, 0.0, 0.0) 
          CALL FPRINTF('based on arsa.f by Brusco, M., '
     1 //'Koehn, H.F.,', 41, 0.0, 0.0) 
          CALL FPRINTF('and Stahl, S. (2007)', 21, 0.0, 0.0) 
          CALL FPRINTF('', 0, 0.0, 0.0) 

C          PRINT *, ''
C          PRINT *, 'COOL =', COOL
C          PRINT *, 'TMIN =', TMIN
C          PRINT *, 'NREPS =', NREPS
C          PRINT *, ''

          CALL FPRINTF('COOL = %5.3f', 12, DBLE(COOL), 0.0) 
          CALL FPRINTF('TMIN = %5.3f', 12, DBLE(TMIN), 0.0) 
          CALL FPRINTF('NREPS= %5.0f', 12, DBLE(NREPS), 0.0) 
          CALL FPRINTF('', 0, 0.0, 0.0) 

      ENDIF

C      INTEGER U(400), S(400), UNSEL, T(100,400), SB(400), Q, GB(400)
C      CHARACTER JNK
C      OPEN(1, FILE = 'DATASET')
C      OPEN(2, FILE = 'OUTPUT')
C      OPEN(3, FILE = 'PERMUT')
C      CALL GETTIM (IHR, IMIN, ISEC, I100)
C      CALL GETDAT (IYR, IMON, IDAY)  
C      TIMEA = DFLOAT(86400*IDAY+3600*IHR+60*IMIN+ISEC)+DFLOAT(I100)/100.  
      RULE = .5
C      READ(1,*) N
C      READ(1,*) JNK
C      READ(1,*) ((A(I,J),J=1,N),I=1,N)
      ICT = 0
      DO I = 1,N-1
        DO J = I+1,N
          ICT = ICT + 1
          D(I,J) = DFLOAT(J-I)
          D(J,I) = D(I,J)
          R1(ICT) = D(I,J)
          R2(ICT) = A(I,J)
        END DO
      END DO
      DO I = 1, ICT-1
        DO J = I+1,ICT
          IF(R1(J).GT.R1(I)) THEN
            RDUM = R1(J)
            R1(J) = R1(I)
            R1(I) = RDUM
          END IF
          IF(R2(J).GT.R2(I)) THEN
            RDUM = R2(J)
            R2(J) = R2(I)
            R2(I) = RDUM
          END IF
        END DO
      END DO
      ASUM = 0.0D0
      DO I = 1,ICT
        ASUM = ASUM + R1(I) * R2(I)
      END DO
      EPS = 1.0D-08
C      NREPS = 20
      DO 999 III = 1,NREPS
        DO I = 1,N
          U(I) = I
          T(III,I) = 0
        END DO
        UNSEL = N
        DO 1 I = 1,N
C          S1 = rand()
          CALL unifrand(S1)
          ISET = S1 * FLOAT(UNSEL) + 1.
          IF(ISET.GT.UNSEL) ISET = UNSEL
          T(III,I) = U(ISET)
C          DO J = ISET,UNSEL   
C	    out of bounds error reported by Rohan Shah (9/13/12)
          DO J = ISET,UNSEL-1
            U(J) = U(J+1)
          END DO
          UNSEL = UNSEL - 1
  1     CONTINUE
 999  CONTINUE
C
      ZMIN = 9.9D+20
      ZAVG = 0.
      ZMAX = 0.
      DO 1000 III = 1,NREPS
        DO I = 1,N
          S(I) = T(III,I)
        END DO
        Z = 0.0D0
        DO I = 1,N-1
          K = S(I)
          DO J = I+1,N
            L = S(J)
            Z = Z + D(I,J) * A(K,L)
          END DO
        END DO
        ZBEST = Z
        TMAX = 0.0D0
        DO LLL = 1,5000
C          S1 = rand()
          CALL unifrand(S1)
          I1 = S1 * FLOAT(N) + 1.
          IF(I1.GT.N) I1 = N
C 199      S1 = rand()
 199      CALL unifrand(S1)
          J1 = S1 * FLOAT(N) + 1.
          IF(J1.GT.N) J1 = N
          IF(I1.EQ.J1) GO TO 199
          IF(I1.GT.J1) THEN
            JDUM = J1
            J1 = I1
            I1 = JDUM
          END IF
          K = S(I1)
          M = S(J1)
          DELTA = 0.0D0
          DO 1250 L1 = 1,N
            IF(I1.EQ.L1.OR.J1.EQ.L1) GO TO 1250
            L=S(L1)
            DELTA=DELTA+(D(L1,I1)-D(L1,J1))*(A(L,M)-A(L,K))
 1250     CONTINUE
          IF(DELTA.LT.0) THEN
            IF(ABS(DELTA).GT.TMAX) TMAX = ABS(DELTA)
          END IF
        END DO
C        COOL = .95
C        TMIN = .0001d0
C        TMAX = Z
        ILOOP = 100*N
        NLOOP = (LOG(TMIN)-LOG(TMAX))/LOG(COOL)
C        WRITE(*,21) TMIN,TMAX,NLOOP
C  21    FORMAT(2F14.5,I6)
        IF (IVERB == 1) THEN
C            WRITE(*,21) NLOOP
            CALL FPRINTF('Steps needed:  %10.0f', 21, DBLE(NLOOP), 0.0) 
        ENDIF
C  21    FORMAT('Steps needed: ',I10)
C        GO TO 889
        TEMP = TMAX
        DO I = 1,N
          SB(I) = S(I)
        END DO
C
        DO 2000 IJK = 1,NLOOP
        IF (IVERB == 1) THEN
C            WRITE(*,22) TEMP
            CALL FPRINTF('Temp = %14.5f', 13, DBLE(TEMP), 0.0) 
        ENDIF
C  22    FORMAT('Temp = ', F14.5)

C   R interrupt
        CALL rchkusr()
C

          DO 2001 KKK = 1,ILOOP
C            S1 = rand()
            CALL unifrand(S1)
            IF(S1.LE.RULE) THEN     ! INTERCHANGE / INSERTION / OR BOTH
C            S1 = rand()
            CALL unifrand(S1)
            I1 = S1 * FLOAT(N) + 1.
            IF(I1.GT.N) I1 = N
C 99         S1 = rand()
 99         CALL unifrand(S1)
            J1 = S1 * FLOAT(N) + 1.
            IF(J1.GT.N) J1 = N
            IF(I1.EQ.J1) GO TO 99
            IF(I1.GT.J1) THEN
              JDUM = J1
              J1 = I1
              I1 = JDUM
            END IF
            K = S(I1)
            M = S(J1)
            DELTA = 0.0D0
            DO 250 L1 = 1,N
              IF(I1.EQ.L1.OR.J1.EQ.L1) GO TO 250
              L=S(L1)
              DELTA=DELTA+(D(L1,I1)-D(L1,J1))*(A(L,M)-A(L,K))
  250       CONTINUE
            IF(DELTA.GT.-EPS) THEN
              Z = Z + DELTA
              S(I1) = M
              S(J1) = K
              IF(Z.GT.ZBEST) THEN
                ZBEST = Z
                DO I = 1,N
                  SB(I) = S(I)
                END DO
              END IF
            ELSE
C              S1 = rand()
              CALL unifrand(S1)
              RCRIT = EXP(DELTA/TEMP)
              IF(S1.LE.RCRIT) THEN
                Z = Z + DELTA
                S(I1) = M
                S(J1) = K
              END IF
            END IF

            ELSE                ! INSERTION

C            S1 = rand()
            CALL unifrand(S1)
            I1 = S1 * FLOAT(N) + 1.      ! OBJECT POSITION IS I1
            IF(I1.GT.N) I1 = N
C 599        S1 = rand()
 599        CALL unifrand(S1)
            J1 = S1 * FLOAT(N) + 1.
            IF(J1.GT.N) J1 = N
            IF(I1.EQ.J1) GO TO 599
            K = S(I1)
            DELTA = 0.0D0
            IF(J1.GT.I1) THEN
              SPAN = DFLOAT(J1-I1)
              DO L = I1+1,J1
                Q = S(L)
                DO I = J1+1,N
                  M = S(I)
                  DELTA = DELTA + A(M,Q)
                END DO
                DO I = 1,I1-1
                  M = S(I)
                  DELTA = DELTA - A(M,Q)
                END DO
              END DO
              DO I = 1,I1-1
                M = S(I)
                DELTA = DELTA + SPAN*A(M,K)
              END DO
              DO I = J1+1,N
                M = S(I)
                DELTA = DELTA - SPAN*A(K,M)
              END DO
              SPAN2 = SPAN+1
              DO I = I1+1,J1
                SPAN2 = SPAN2-2
                M = S(I)
                DELTA = DELTA + SPAN2*A(K,M)
              END DO
            ELSE
              SPAN = DFLOAT(I1-J1)
              DO L = J1,I1-1
                Q = S(L)
                DO I = I1+1,N
                  M = S(I)
                  DELTA = DELTA - A(M,Q)
                END DO
                DO I = 1,J1-1
                  M = S(I)
                  DELTA = DELTA + A(M,Q)
                END DO
              END DO
              DO I = 1,J1-1
                M = S(I)
                DELTA = DELTA - SPAN*A(M,K)
              END DO
              DO I = I1+1,N
                M = S(I)
                DELTA = DELTA + SPAN*A(K,M)
              END DO
              SPAN2 = SPAN+1
              DO I = J1,I1-1
                SPAN2 = SPAN2-2
                M = S(I)
                DELTA = DELTA - SPAN2*A(K,M)
              END DO
            END IF  
            IF(DELTA.GT.-EPS) THEN
              Z = Z + DELTA
              IF(J1.GT.I1) THEN
                DO L = I1,J1-1
                  S(L)=S(L+1)
                END DO
                S(J1) = K
              ELSE
                DO L = I1,J1+1,-1
                  S(L)=S(L-1)
                END DO
                S(J1) = K
              END IF

              IF(Z.GT.ZBEST) THEN
                ZBEST = Z
                DO I = 1,N
                  SB(I) = S(I)
                END DO
              END IF
            ELSE
C              S1 = rand()
              CALL unifrand(S1)
              RCRIT = EXP(DELTA/TEMP)
              IF(S1.LE.RCRIT) THEN
                Z = Z + DELTA
                IF(J1.GT.I1) THEN
                  DO L = I1,J1-1
                    S(L)=S(L+1)
                  END DO
                  S(J1) = K
                ELSE
                  DO L = I1,J1+1,-1
                    S(L)=S(L-1)
                  END DO
                  S(J1) = K
                END IF
              END IF
            END IF
C
            END IF
 2001     CONTINUE
          TEMP = TEMP*COOL
 2000   CONTINUE
C
C        DO I = 1,N
C          WRITE(3,*) SB(I)
C        END DO
C        SOLS(III) = ZBEST
C        RMED(III) = ZBEST
C        ZSUM = ZSUM + ZBEST
        IF(ZBEST.LT.ZMIN) ZMIN = ZBEST
        IF(ZBEST.GT.ZMAX) THEN
          ZMAX = ZBEST
          DO I = 1,N
            IPERM(I) = SB(I)
          END DO
        END IF
C        WRITE(*,77) iii,zbest
 1000 CONTINUE
C      CALL GETTIM (IHR, IMIN, ISEC, I100)  
C      CALL GETDAT (IYR, IMON, IDAY)  
C      TIMEB = DFLOAT(86400*IDAY+3600*IHR+60*IMIN+ISEC)+DFLOAT(I100)/100.  
C      TIMTOT = TIMEB - TIMEA
C      ZSUM = ZSUM/DFLOAT(NREPS)
C      TIMTOT = TIMTOT / DFLOAT(NREPS)
C      IBEST = 0
C      DO I = 1,NREPS
C        IF(SOLS(I).GT.ZMAX-1) IBEST = IBEST + 1
C      END DO
C      ATTR = 100.*(DFLOAT(IBEST) / DFLOAT(NREPS))
C      DO I = 1,NREPS-1
C        DO J = I + 1,NREPS
C          IF(RMED(J).LT.RMED(I)) THEN
C            RST = RMED(J)
C            RMED(J) = RMED(I)
C            RMED(I) = RST
C          END IF
C        END DO
C      END DO
C      NR = NREPS/2
C      RMN = (RMED(NR)+RMED(NR+1))/2.0
C      WRITE(3,*) N
C      DO I = 1,N
C        WRITE(3,*) I,GB(I)
C      END DO
C      WRITE(2,201)
C      WRITE(*,201)
C      WRITE(2,202) ZMIN, ZSUM, RMN, ZMAX,ATTR,TIMTOT
C      WRITE(*,202) ZMIN, ZSUM, RMN, ZMAX,ATTR,TIMTOT
C      WRITE(2,203) ZMAX,ZMAX/ASUM,ATTR,TIMTOT/20.
C      WRITE(*,203) ZMAX,ZMAX/ASUM,ATTR,TIMTOT/20.

C    Return R RNG
      CALL Putrngstate()

      RETURN
   
C   77 format(I5,f15.4)
C  200 FORMAT(I3,2F20.6)
C  201 FORMAT(9X,'MIN',11X,'AVG',11X,'MED',10X,'MAX',5X,'ATT',3X,'TIME')
C  202 FORMAT(4F18.2,F5.0,F7.2)
C  203 FORMAT(F20.4,F11.8,F5.0,F7.2)
  889 CONTINUE
      END

C   ANTI-ROBINSON SERIATION
C   simulated annealing algorithm - provides an initial permutation
C   by Brusco, M., Koehn, H.F., and Stahl, S.
C   R Interface by Michael Hahsler

C      PROGRAM SANNEAL
      SUBROUTINE arsa(N, A, COOL, TMIN, NREPS, IPERM, D, U,
     1 S, T, SB, ZMAX, RULE, TRYMULT, IVERB)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N,N)
      DIMENSION IPERM(N)
      DOUBLE PRECISION D(N,N)
      REAL S1, RCRIT
      INTEGER U(N), S(N), UNSEL, T(NREPS,N), SB(N), Q, NREPS

      EPS = 1.0D-08

C   Defaults
C      RULE = .5
C        COOL = .95
C        TMIN = .0001d0

C     Initialize R RNG
      CALL getrngstate()

      DO I = 1,N-1
        DO J = I+1,N
          D(I,J) = DFLOAT(J-I)
          D(J,I) = D(I,J)
        END DO
      END DO

      DO 999 III = 1,NREPS
        DO I = 1,N
          U(I) = I
          T(III,I) = 0
        END DO
        UNSEL = N
        DO 1 I = 1,N
C          S1 = rand()
          CALL unifrand(S1)
          ISET = INT(S1 * FLOAT(UNSEL) + 1.)
          IF(ISET.GT.UNSEL) ISET = UNSEL
          T(III,I) = U(ISET)
C          DO J = ISET,UNSEL
C          out of bounds error reported by Rohan Shah (9/13/12)
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
C        DO LLL = 1,5000
C   Find initial TMAX using N*10 tries
        DO LLL = 1,N*10
C          S1 = rand()
          CALL unifrand(S1)
          I1 = INT(S1 * FLOAT(N) + 1.)
          IF(I1.GT.N) I1 = N
C 199      S1 = rand()
 199      CALL unifrand(S1)
          J1 = INT(S1 * FLOAT(N) + 1.)
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
C        TMAX = Z
        ILOOP = INT(TRYMULT*N)
        NLOOP = INT((LOG(TMIN)-LOG(TMAX))/LOG(COOL))
        IF (IVERB == 1) THEN
            CALL intpr('Steps needed', -1, NLOOP, 1)
            CALL intpr('Temp', -1, NLOOP, 0)
        ENDIF
        TEMP = TMAX
        DO I = 1,N
          SB(I) = S(I)
        END DO
C
        DO 2000 IJK = 1,NLOOP
        IF (IVERB == 1) THEN
            CALL dblepr('', -1, DBLE(TEMP), 1)
        ENDIF

C   R interrupt
        CALL rchkusr()
C

          DO 2001 KKK = 1,ILOOP
C            S1 = rand()
            CALL unifrand(S1)
            IF(S1.LE.RULE) THEN     ! INTERCHANGE / INSERTION / OR BOTH
C            S1 = rand()
            CALL unifrand(S1)
            I1 = INT(S1 * FLOAT(N) + 1.)
            IF(I1.GT.N) I1 = N
C 99         S1 = rand()
 99         CALL unifrand(S1)
            J1 = INT(S1 * FLOAT(N) + 1.)
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
              RCRIT = REAL(EXP(DELTA/TEMP))
              IF(S1.LE.RCRIT) THEN
                Z = Z + DELTA
                S(I1) = M
                S(J1) = K
              END IF
            END IF

            ELSE                ! INSERTION

C            S1 = rand()
            CALL unifrand(S1)
            I1 = INT(S1 * FLOAT(N) + 1.)      ! OBJECT POSITION IS I1
            IF(I1.GT.N) I1 = N
C 599        S1 = rand()
 599        CALL unifrand(S1)
            J1 = INT(S1 * FLOAT(N) + 1.)
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
              RCRIT = REAL(EXP(DELTA/TEMP))
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

        IF(ZBEST.LT.ZMIN) ZMIN = ZBEST
        IF(ZBEST.GT.ZMAX) THEN
          ZMAX = ZBEST
          DO I = 1,N
            IPERM(I) = SB(I)
          END DO
        END IF
        IF (IVERB == 1) THEN
          CALL intpr('Rep', -1, III, 1)
          CALL dblepr('ZMAX', -1, DBLE(ZMAX), 1)
        END IF
 1000 CONTINUE

C    Return R RNG
      CALL Putrngstate()

      RETURN
      END

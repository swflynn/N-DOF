      BLOCK DATA BDSOBL
      INTEGER POLY(2:40),VINIT(2:40,8)
      COMMON  /SOBDAT/ POLY,VINIT
      SAVE    /SOBDAT/
      DATA POLY    /3,7,11,13,19,25,37,59,47,
     +              61,55,41,67,97,91,109,103,115,131,
     +              193,137,145,143,241,157,185,167,229,171,
     +              213,191,253,203,211,239,247,285,369,299/
      DATA (VINIT(I,1),I=2,40)  /39*1/
      DATA (VINIT(I,2),I=3,40)  /1,3,1,3,1,3,3,1,
     +                           3,1,3,1,3,1,1,3,1,3,
     +                           1,3,1,3,3,1,3,1,3,1,
     +                           3,1,1,3,1,3,1,3,1,3/
      DATA (VINIT(I,3),I=4,40)  /7,5,1,3,3,7,5,
     +                           5,7,7,1,3,3,7,5,1,1,
     +                           5,3,3,1,7,5,1,3,3,7,
     +                           5,1,1,5,7,7,5,1,3,3/
      DATA (VINIT(I,4),I=6,40)  /1,7,9,13,11,
     +                           1,3,7,9,5,13,13,11,3,15,
     +                           5,3,15,7,9,13,9,1,11,7,
     +                           5,15,1,15,11,5,3,1,7,9/
      DATA (VINIT(I,5),I=8,40)  /9,3,27,
     +                           15,29,21,23,19,11,25,7,13,17,
     +                           1,25,29,3,31,11,5,23,27,19,
     +                           21,5,1,17,13,7,15,9,31,9/
      DATA (VINIT(I,6),I=14,40) /37,33,7,5,11,39,63,
     +                           27,17,15,23,29,3,21,13,31,25,
     +                           9,49,33,19,29,11,19,27,15,25/
      DATA (VINIT(I,7),I=20,40) /13,
     +                           33,115,41,79,17,29,119,75,73,105,
     +                           7,59,65,21,3,113,61,89,45,107/
      DATA (VINIT(I,8),I=38,40) /7,23,39/
      END
      
      SUBROUTINE INSSOBL (FLAG, DIMEN, ATMOST, TAUS, QUASI,
     *                    MAX,IFLAG)

      INTEGER  POLY(2:40),VINIT(2:40,8),MAX
      INTEGER  V(40,31),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  DIMEN,ATMOST,I,J,K,P,M,NEWV,EXOR,L
      INTEGER  TAU(13),TAUS,SHIFT(40),LSM(40,31),SV(40,31)
      INTEGER  OUTPUT(22),DIGIT(60),TEMP1,TEMP2,TEMP3,TEMP4
      INTEGER  USM(31,31),USHIFT(31),IFLAG,TV(40,31,31),MAXX
      REAL*8   RECIPD,QUASI(40),LL
      LOGICAL  FLAG(2),INCLUD(8)
      EXTERNAL EXOR
      COMMON   /SOBDAT/ POLY,VINIT
      COMMON   /SOBOL/  S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /SOBDAT/,/SOBOL/
      DATA TAU /0,0,1,3,5,8,11,15,19,23,27,31,35/

      S = DIMEN
      FLAG(1) = (S .GE. 1 .AND. S .LE. 40)
      FLAG(2) = (ATMOST .GT.0 .AND. ATMOST .LT. 2**30)
      IF (.NOT. (FLAG(1) .AND. FLAG(2))) RETURN
      IF (S .LE. 13) THEN
        TAUS = TAU(S)
      ELSE
        TAUS = -1
      ENDIF

      I = ATMOST
      MAXCOL = 0
   10 MAXCOL = MAXCOL + 1
      I = I / 2
      IF (I .GT. 0) GOTO 10
      

      DO 20 I = 1, MAXCOL
   20 V(1,I) = 1

      DO 100 I = 2, S

        J = POLY(I)
        M = 0
   30   J = J / 2
        IF (J .GT. 0) THEN
          M = M + 1
          GOTO 30
        ENDIF

        J = POLY(I)
        DO 40 K = M, 1, -1
          INCLUD(K) = (MOD(J,2) .EQ. 1)
          J = J / 2
   40   CONTINUE

        DO 50 J = 1, M
          V(I,J) = VINIT(I, J)
   50   CONTINUE

        DO 70 J = M+1, MAXCOL
          NEWV = V(I, J-M)
          L = 1
          DO 60 K = 1, M
            L = 2 * L
            IF (INCLUD(K)) NEWV = EXOR(NEWV, L * V(I, J-K))

   60     CONTINUE
          V(I,J) = NEWV
          
   70   CONTINUE
      
  100 CONTINUE

      L = 1
      DO 120 J = MAXCOL-1, 1, -1
        L = 2 * L
       DO 110 I = 1, S
           V(I,J) = V(I,J) * L
  110   CONTINUE
  120 CONTINUE

      IF (IFLAG .EQ. 0) THEN
         DO 130 I = 1,S
            DO 140 J= 1,MAXCOL 
                  SV(I,J) = V(I,J)
 140          CONTINUE
              SHIFT(I) = 0
 130      CONTINUE
          LL= 2.0**(MAXCOL)
      ELSE             
        IF ((IFLAG .EQ. 1) .OR. (IFLAG .EQ. 3)) THEN
         CALL GENSCRML(MAX,LSM,SHIFT)
         DO  230 I = 1,S
           DO 240 J = 1,MAXCOL
             L = 1
             TEMP2 = 0
             DO 250 P = MAX,1,-1
                TEMP1 = 0
                DO 260 K = 1,MAXCOL
                   TEMP1 = TEMP1+ 
     *                   (IBITS(LSM(I,P),K-1,1)*IBITS(V(I,J),K-1,1))                                         
 260            CONTINUE
                TEMP1 = MOD(TEMP1,2)
                TEMP2 = TEMP2+TEMP1*L   
                L = 2 * L
 250          CONTINUE
              SV(I,J) = TEMP2
 240      CONTINUE
 230     CONTINUE
            LL= 2.0**(MAX)
       ENDIF
       IF ((IFLAG .EQ. 2) .OR. (IFLAG .EQ. 3)) THEN
        CALL GENSCRMU(USM,USHIFT) 
         IF (IFLAG .EQ. 2) THEN
           MAXX = MAXCOL
         ELSE
           MAXX = MAX
         ENDIF    
        DO 330 I = 1,S
         DO 340 J = 1,MAXCOL
          P = MAXX
          DO 350 K = 1,MAXX
             IF (IFLAG .EQ. 2) THEN
              TV(I,P,J) = IBITS(V(I,J),K-1,1)
             ELSE
               TV(I,P,J) = IBITS(SV(I,J),K-1,1) 
             ENDIF 
              P = P-1
 350        CONTINUE
 340       CONTINUE
        
          DO 341 PP = 1,MAXCOL 
            TEMP2 = 0 
            TEMP4 = 0
            L = 1
            DO 360 J = MAXX,1,-1
                TEMP1 = 0
                TEMP3 = 0
             DO 370 P = 1,MAXCOL
                TEMP1 = TEMP1 + TV(I,J,P)*USM(P,PP)
                IF (PP .EQ. 1) THEN
                  TEMP3 = TEMP3 + TV(I,J,P)*USHIFT(P)
                ENDIF 
370           CONTINUE
              TEMP1 = MOD(TEMP1,2)
              TEMP2 = TEMP2 + TEMP1*L
              IF (PP .EQ. 1) THEN 
                TEMP3  = MOD(TEMP3,2)
                TEMP4 = TEMP4 + TEMP3*L
              ENDIF  
              L = 2*L
 360         CONTINUE
              SV(I,PP) = TEMP2
              IF (PP .EQ. 1) THEN
                 IF (IFLAG .EQ. 3) THEN
                   SHIFT(I) = EXOR(TEMP4, SHIFT(I))           
                 ELSE
                   SHIFT(I) = TEMP4
                 ENDIF  
              ENDIF
 341        CONTINUE 
 330       CONTINUE
             LL = 2.0**(MAXX)
       ENDIF
      ENDIF 

      RECIPD = 1.0 /LL

      COUNT = 0
      DO 400 I = 1, S
        LASTQ(I) = SHIFT(I)
        QUASI(I) = LASTQ(I)*RECIPD
 400   CONTINUE
      RETURN
      END

      SUBROUTINE GENSCRML(MAX,LSM,SHIFT)

      INTEGER LSM(40,31),MAXCOL,P,I,J
      INTEGER SHIFT(40),MAX,S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /SOBOL/ S,MAXCOL
      SAVE /SOBOL/
      
      DO 10 P = 1,S
               SHIFT(P) = 0
               L = 1
         DO 20 I = MAX,1,-1
               LSM(P,I) = 0
               STEMP =  MOD((int(UNI()*1000.0)),2)
               SHIFT(P) = SHIFT(P)+STEMP*L
               L = 2 * L
               LL = 1
            DO 30 J = MAXCOL,1,-1
               IF (J .EQ. I) THEN
                TEMP = 1
               ELSE IF (J .LT. I)  THEN 
                TEMP = MOD((int(UNI()*1000.0)),2)
               ELSE
                TEMP = 0
               ENDIF
               LSM(P,I) = LSM(P,I)+TEMP*LL
               LL = 2 * LL            
 30           CONTINUE
 20        CONTINUE 
 10      CONTINUE
      RETURN
      END   
      
      SUBROUTINE GENSCRMU(USM,USHIFT)

      INTEGER USM(31,31),MAXCOL,P,I,J
      INTEGER USHIFT(31),S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /SOBOL/ S,MAXCOL
      SAVE /SOBOL/
      
          DO 20 I = 1,MAXCOL
             STEMP =  MOD((int(UNI()*1000.0)),2)
             USHIFT(I) = STEMP               
             DO 30 J = 1,MAXCOL
               IF (J .EQ. I) THEN
                 TEMP = 1
               ELSE IF (J .GT. I)  THEN 
                 TEMP = MOD((int(UNI()*1000.0)),2)
               ELSE
                TEMP = 0
               ENDIF
               USM(I,J) = TEMP        
 30           CONTINUE
 20        CONTINUE 
      RETURN
      END   
      
      DOUBLE PRECISION FUNCTION UNI()
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY
      PARAMETER ( TWOM24 = 1D0/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS / 
     & 0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,
     & 0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043,
     & 0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127,
     & 0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      
      
      UNI = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNI .LT. 0 ) THEN 
         UNI = UNI + 1
         CARRY = TWOM24
      ELSE 
         CARRY = 0
      ENDIF
      SEEDS(I) = UNI
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )
      RETURN
      END      

      SUBROUTINE GOSSOBL (QUASI)

      INTEGER  SV(40,31),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  I,L,EXOR
      REAL*8     RECIPD,QUASI(40)
      COMMON   /SOBOL/ S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /SOBOL/

      L = 0
      I = COUNT
    1 L = L + 1
      IF (MOD(I,2) .EQ. 1) THEN
        I = I / 2
        GOTO 1
      ENDIF
      IF (L .GT. MAXCOL) STOP ' TOO MANY CALLS ON GOSOBL'
      DO 2 I = 1, S
        LASTQ(I) = EXOR (LASTQ(I), SV(I,L))
        QUASI(I) = LASTQ(I) * RECIPD
    2 CONTINUE
      COUNT = COUNT + 1
      RETURN
      END
      
      INTEGER FUNCTION EXOR (IIN, JIN)
      INTEGER I,J,K,L,IIN,JIN

      I = IIN
      J = JIN
      K = 0
      L = 1
    1 IF (I .EQ. J) THEN
        EXOR = K
        RETURN
      ENDIF

      IF (MOD(I,2) .NE. MOD(J,2)) K = K + L
      I = I / 2
      J = J / 2
      L = 2 * L
      GOTO 1
      END

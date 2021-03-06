      SUBROUTINE JSTNUM(STR,SIGN,IRET)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    JSTNUM
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE REMOVES ALL LEADING BLANKS FROM A CHARACTER
C   STRING CONTAINING AN ENCODED INTEGER VALUE.  IF THE VALUE HAS A
C   LEADING SIGN CHARACTER ('+' OR '-'), THEN THIS CHARACTER IS ALSO
C   REMOVED AND IS RETURNED SEPARATELY WITHIN SIGN.  IF THE RESULTANT
C   STRING CONTAINS ANY NON-NUMERIC CHARACTERS, THAN AN APPROPRIATE
C   CALL IS MADE TO TO BUFR ARCHIVE LIBRARY SUBROUTINE BORT.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR (ENTRY POINT IN JSTIFY)
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT" (IN PARENT ROUTINE JSTIFY)
C 2002-05-14  J. WOOLLEN -- CHANGED FROM AN ENTRY POINT TO INCREASE
C                           PORTABILITY TO OTHER PLATFORMS (JSTIFY WAS
C                           THEN REMOVED BECAUSE IT WAS JUST A DUMMY
C                           ROUTINE WITH ENTRIES)
C 2003-11-04  J. ATOR    -- ADDED DOCUMENTATION
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED HISTORY
C                           DOCUMENTATION; OUTPUTS MORE COMPLETE
C                           DIAGNOSTIC INFO WHEN ROUTINE TERMINATES
C                           ABNORMALLY OR UNUSUAL THINGS HAPPEN
C DART $Id$
C
C USAGE:    CALL JSTNUM (STR, SIGN, IRET)
C   INPUT ARGUMENT LIST:
C     STR      - CHARACTER*(*): STRING CONTAINING ENCODED INTEGER VALUE
C
C   OUTPUT ARGUMENT LIST:
C     STR      - CHARACTER*(*): COPY OF INPUT STR WITH LEADING BLANKS
C                AND SIGN CHARACTER REMOVED
C     SIGN     - CHARACTER*1: SIGN OF ENCODED INTEGER VALUE:
C                     '+' = positive value
C                     '-' = negative value
C     IRET     - INTEGER: RETURN CODE:
C                       0 = normal return
C                      -1 = encoded value within STR was not an integer
C
C   OUTPUT FILES:
C     UNIT 06  - STANDARD OUTPUT PRINT
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     STRNUM
C    THIS ROUTINE IS CALLED BY: ELEMDX
C                               Normally not called by any application
C                               programs but it could be.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      CHARACTER*(*) STR
      CHARACTER*1  SIGN

      COMMON /QUIET / IPRT

      IRET = 0

      IF(STR.EQ.' ') GOTO 900

      LSTR = LEN(STR)
2     IF(STR(1:1).EQ.' ') THEN
         STR  = STR(2:LSTR)
         GOTO 2
      ENDIF
      IF(STR(1:1).EQ.'+') THEN
         STR  = STR(2:LSTR)
         SIGN = '+'
      ELSEIF(STR(1:1).EQ.'-') THEN
         STR  = STR(2:LSTR)
         SIGN = '-'
      ELSE
         SIGN = '+'
      ENDIF

      CALL STRNUM(STR,NUM)
      IF(NUM.LT.0) THEN
         IF(IPRT.GE.0) THEN
      PRINT*
      PRINT*,'+++++++++++++++++++++++WARNING+++++++++++++++++++++++++'
         PRINT*, 'BUFRLIB: JSTNUM: ENCODED VALUE WITHIN RESULTANT ',
     .    'CHARACTER STRING (',STR,') IS NOT AN INTEGER - RETURN WITH ',
     .    'IRET = -1'
      PRINT*,'+++++++++++++++++++++++WARNING+++++++++++++++++++++++++'
      PRINT*
         ENDIF
         IRET = -1
      ENDIF

C  EXITS
C  -----

      RETURN
900   CALL BORT('BUFRLIB: JSTNUM - INPUT BLANK CHARACTER STRING NOT '//
     . 'ALLOWED')
      END

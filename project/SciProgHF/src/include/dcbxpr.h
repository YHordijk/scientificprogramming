!
!     PROPERTIES:
!
!       IPRPTYP  - indicates kind of property operator
!         0. FOCKMAT       * Fock-matrix
!         1. P             * scalar operator
!         2. [alpha_x]p    * x-component of alpha times scalar operator
!         3. [alpha_y]p    * y-component of alpha times scalar operator
!         4. [alpha_z]p    * z-component of alpha times scalar operator
!         5. [Alpha x P]_x * vector product of alpha and vector operator,
!                            x-component    
!         6. [Alpha x P]_y * vector product of alpha and vector operator,
!                            y-component    
!         7. [Alpha x P]_z * vector product of alpha and vector operator,
!                            z-component    
!         8. A.P           * dot-product of alpha and vector operator
!         9. gamma5 p      * gamma5 times scalar operator
!        10. [sigma_x]p    * x-component of sigma times scalar operator 
!        11. [sigma_y]p    * y-component of sigma times scalar operator
!        12. [sigma_z]p    * z-component of sigma times scalar operator
!        13. [betasig_x]p  * x-component of beta sigma times scalar operator
!        14. [betasig_y]p  * y-component of beta sigma times scalar operator
!        15. [betasig_z]p  * z-component of beta sigma times scalar operator
!        16. [betaalp_x]p  * x-component of beta alpha times scalar operator
!        17. [betaalp_y]p  * y-component of beta alpha times scalar operator
!        18. [betaalp_z]p  * z-component of beta alpha times scalar operator 
!        20. S.P           * dot-product of sigma and vector operator
!
!       IPRPLBL - pointers to property labels
!       IPRPSYM - boson symmetry of operator  (1 - 8)
!       IPRPTIM - symmetry under time reversal (+1,-1)
!       PRPNAM  - user specified name of operator
!
      INTEGER MAXPRPS
      PARAMETER (MAXPRPS = 5328)

      CHARACTER*16 PRPNAM(MAXPRPS)
      REAL*8 FACPRP(3, MAXPRPS)
      INTEGER NPRPS, IPRPTYP(MAXPRPS), IPRPLBL(3, MAXPRPS),             &
     &        IPRPSYM(MAXPRPS), IPRPTIM(MAXPRPS)
      COMMON /XCBPRP/ FACPRP, NPRPS, IPRPTYP, IPRPLBL, IPRPSYM, IPRPTIM

      COMMON /XCBPRPC/ PRPNAM

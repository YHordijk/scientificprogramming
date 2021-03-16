#if defined (SYS_CRAY) || defined (SYS_T3D)
!DIR$ IVDEP
#endif
#if defined (SYS_SX)
!VDIR NODEP
#endif
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (SYS_SX)
!vectorization note: ignore vector dependence
#endif

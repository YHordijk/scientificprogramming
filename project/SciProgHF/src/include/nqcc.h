      real(8) ::    efgel     (3, 3, mxcent),                           &
     &              efgnuc    (3, 3, mxcent),                           &
     &              efgcorr   (3, 3, mxcent),                           &
     &              efgtot    (3, 3, mxcent),                           &
     &              efg_method(3, 3, mxcent)

      common /nqcc/ efgel,                                              &
     &              efgnuc,                                             &
     &              efgcorr,                                            &
     &              efgtot,                                             &
     &              efg_method

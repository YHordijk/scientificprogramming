*Comdeck Regge $Revision$
      Integer facul,prim,nprim,iwork 
      parameter (nprim=11,mxLinRE=36)
cbs   nprim is the number of prime-numbers
      dimension facul(nprim,0:mxLinRE),prim(nprim),
     *iwork(nprim),ihigh(0:mxLinRE)
      data prim /2,3,5,7,11,13,17,19,23,29,31/        !prime numbers 
c     
c     decompose facultatives into powers of prime numbers  
c    
      Data facul / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     >             0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0,
     >             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     >             1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     >             3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     >             3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     >             4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     >             4, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     >             7, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     >             7, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     >             8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0,
     >             8, 4, 2, 1, 1, 0, 0, 0, 0, 0, 0,
     >            10, 5, 2, 1, 1, 0, 0, 0, 0, 0, 0,
     >            10, 5, 2, 1, 1, 1, 0, 0, 0, 0, 0,
     >            11, 5, 2, 2, 1, 1, 0, 0, 0, 0, 0, 
     >            11, 6, 3, 2, 1, 1, 0, 0, 0, 0, 0, 
     >            15, 6, 3, 2, 1, 1, 0, 0, 0, 0, 0,
     >            15, 6, 3, 2, 1, 1, 1, 0, 0, 0, 0,
     >            16, 8, 3, 2, 1, 1, 1, 0, 0, 0, 0, 
     >            16, 8, 3, 2, 1, 1, 1, 1, 0, 0, 0,
     >            18, 8, 4, 2, 1, 1, 1, 1, 0, 0, 0, 
     >            18, 9, 4, 3, 1, 1, 1, 1, 0, 0, 0, 
     >            19, 9, 4, 3, 2, 1, 1, 1, 0, 0, 0, 
     >            19, 9, 4, 3, 2, 1, 1, 1, 1, 0, 0,
     >            22,10, 4, 3, 2, 1, 1, 1, 1, 0, 0, 
     >            22,10, 6, 3, 2, 1, 1, 1, 1, 0, 0, 
     >            23,10, 6, 3, 2, 2, 1, 1, 1, 0, 0, 
     >            23,13, 6, 3, 2, 2, 1, 1, 1, 0, 0, 
     >            25,13, 6, 4, 2, 2, 1, 1, 1, 0, 0,
     >            25,13, 6, 4, 2, 2, 1, 1, 1, 1, 0,
     >            26,14, 7, 4, 2, 2, 1, 1, 1, 1, 0, 
     >            26,14, 7, 4, 2, 2, 1, 1, 1, 1, 1, 
     >            31,14, 7, 4, 2, 2, 1, 1, 1, 1, 1, 
     >            31,15, 7, 4, 3, 2, 1, 1, 1, 1, 1, 
     >            32,15, 7, 4, 3, 3, 1, 1, 1, 1, 1, 
     >            32,15, 8, 5, 3, 3, 1, 1, 1, 1, 1, 
     >            34,17, 8, 5, 3, 3, 1, 1, 1, 1, 1/ 
c
      data ihigh /0,0,1,2,2,3,3,4,4,4,4,5,5,
     >            6,6,6,6,7,7,8,8,8,8,9,9,9,
     >            9,9,9,10,10,11,11,11,11,11,11/

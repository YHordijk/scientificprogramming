 module symmetry_offset
 implicit none 

  public auto_symmetry_offset
  public auto_symmetry_offset_triangular
  public alloc_array 
  public dealloc_array 

 private
    public Offset
    type Offset
      integer,allocatable :: oneDirac(:),oneNonDirac(:)
      integer,allocatable :: twoDirac(:,:),twoNonDirac(:,:)
    endtype Offset

  contains

  subroutine auto_symmetry_offset(a,row,column,shift_for_boson_row,shift_for_boson_column)

#include "symm.inc"
! In this subroutine my target is to generate the offset arrays for symmetry.  
   integer,intent(in)  :: row(nrep),column(nrep)
   integer             :: irep, jrep, ijrep
   logical,intent(in)  :: shift_for_boson_row, shift_for_boson_column! a logical variable indicating whether the shift for bosonic irrreps is necessary 
   type(offset),intent(inout) :: a
   integer :: shift1,shift2

! 1> case 1: Offset for (2,2) type integrals.
!            here we have dealt with four different type of offsets 
!            oneDirac ->  bra(k)*ket(l) = ket(kl) !multb(*,*,2)
!                      or,ket(k)*bra(l) = bra(kl)
!            twoDirac ->  bra(k)*ket(l) = ket(k,l)!multb(*,*,2) 
!                      or,ket(k)*bra(l) = bra(k,l)

!            oneNonDirac ->  bra(k)*bra(l) = bra(kl) !multb(*,*,1)
!                      or,ket(k)*ket(l) = ket(kl)
!            twoNonDirac ->  bra(k)*bra(l) = bra(k,l)!multb(*,*,1) 
!                      or,ket(k)*ket(l) = ket(k,l)

!  If the irreps are bosonic we have used a shift of nrep (the total number of irreps) 

   shift1 = 0
   shift2 = 0
   if(shift_for_boson_row)    shift1=nrep 
   if(shift_for_boson_column) shift2=nrep 

   a%oneDirac(1:nrep) = 0
   a%oneNonDirac(1:nrep) = 0
   a%twoDirac(1:nrep,1:nrep) = 0
   a%twoNonDirac(1:nrep,1:nrep) = 0

    do jrep = 1,nrep
      do irep = 1,nrep
        ijrep = multb(irep+shift1,jrep+shift2,2)
        a%twoDirac(irep,jrep) = a%oneDirac(ijrep) 
        a%oneDirac(ijrep) = a%oneDirac(ijrep)+row(irep)*column(jrep) 

        ijrep = multb(irep+shift1,jrep+shift2,1)
        a%twoNonDirac(irep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(irep)*column(jrep) 
      enddo 
    enddo 

   end subroutine

!=============================================================================
  subroutine auto_symmetry_offset_triangular(a,row,column)

#include "symm.inc"
! In this subroutine my target is to generate the offset arrays for symmetry.  
   integer,intent(in)  :: row(nrep),column(nrep)
   integer             :: irep, jrep, ijrep
   type(Offset),intent(inout) :: a


   a%oneNonDirac(1:nrep) = 0
   a%twoNonDirac(1:nrep,1:nrep) = 0

    do jrep = 1,nrep
        ijrep = multb(jrep,jrep,1)
        a%twoNonDirac(jrep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(jrep)*(column(jrep)-1)/2 

      do irep = jrep+1, nrep
        ijrep = multb(irep,jrep,1)
        a%twoNonDirac(irep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(irep)*column(jrep) 
      enddo 
    enddo 

  end subroutine

  subroutine alloc_array(c,array_dim,d,e,f,g)

   type(Offset),intent(inout) :: c
   type(Offset),intent(inout),optional :: d
   type(Offset),intent(inout),optional :: e
   type(Offset),intent(inout),optional :: f
   type(Offset),intent(inout),optional :: g
   integer, intent(in) :: array_dim

   allocate(c%oneDirac(array_dim)) 
   allocate(c%oneNonDirac(array_dim)) 
   allocate(c%twoDirac(array_dim,array_dim)) 
   allocate(c%twoNonDirac(array_dim,array_dim)) 


   if (present(d)) then

   allocate(d%oneDirac(array_dim)) 
   allocate(d%oneNonDirac(array_dim)) 
   allocate(d%twoDirac(array_dim,array_dim)) 
   allocate(d%twoNonDirac(array_dim,array_dim)) 

   endif

  if (present(e)) then

   allocate(e%oneDirac(array_dim)) 
   allocate(e%oneNonDirac(array_dim)) 
   allocate(e%twoDirac(array_dim,array_dim)) 
   allocate(e%twoNonDirac(array_dim,array_dim)) 

   endif
  if (present(f)) then

   allocate(f%oneDirac(array_dim)) 
   allocate(f%oneNonDirac(array_dim)) 
   allocate(f%twoDirac(array_dim,array_dim)) 
   allocate(f%twoNonDirac(array_dim,array_dim)) 

   endif
   if (present(g)) then

   allocate(g%oneDirac(array_dim)) 
   allocate(g%oneNonDirac(array_dim)) 
   allocate(g%twoDirac(array_dim,array_dim)) 
   allocate(g%twoNonDirac(array_dim,array_dim)) 

   endif

  end subroutine

  subroutine dealloc_array(c,d,e,f,g)

   type(Offset),intent(inout) :: c
   type(Offset),intent(inout),optional :: d
   type(Offset),intent(inout),optional :: e
   type(Offset),intent(inout),optional :: f
   type(Offset),intent(inout),optional :: g

  if(allocated(c%oneDirac))    deallocate(c%oneDirac) 
  if(allocated(c%oneNonDirac)) deallocate(c%oneNonDirac) 
  if(allocated(c%twoDirac))    deallocate(c%twoDirac) 
  if(allocated(c%twoNonDirac)) deallocate(c%twoNonDirac) 


   if(present(d)) then

   if(allocated(d%oneDirac))   deallocate(d%oneDirac) 
   if(allocated(d%oneNonDirac))deallocate(d%oneNonDirac) 
   if(allocated(d%twoDirac))   deallocate(d%twoDirac) 
   if(allocated(d%twoNonDirac))deallocate(d%twoNonDirac) 

   endif

   if(present(e)) then

   if(allocated(e%oneDirac))   deallocate(e%oneDirac) 
   if(allocated(e%oneNonDirac))deallocate(e%oneNonDirac) 
   if(allocated(e%twoDirac))   deallocate(e%twoDirac) 
   if(allocated(e%twoNonDirac))deallocate(e%twoNonDirac) 

   endif

   if(present(f)) then

   if(allocated(f%oneDirac))   deallocate(f%oneDirac) 
   if(allocated(f%oneNonDirac))deallocate(f%oneNonDirac) 
   if(allocated(f%twoDirac))   deallocate(f%twoDirac) 
   if(allocated(f%twoNonDirac))deallocate(f%twoNonDirac) 

   endif

   if(present(g)) then

   if(allocated(g%oneDirac))   deallocate(g%oneDirac) 
   if(allocated(g%oneNonDirac))deallocate(g%oneNonDirac) 
   if(allocated(g%twoDirac))   deallocate(g%twoDirac) 
   if(allocated(g%twoNonDirac))deallocate(g%twoNonDirac) 

   endif

  end subroutine

  end module

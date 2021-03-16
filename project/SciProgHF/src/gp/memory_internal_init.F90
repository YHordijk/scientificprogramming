!*
!*
!* Copyright (c) 2008-2013 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*
!
!

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module allocator_internal_init_C
      
         use allocator_parameters

         implicit none
         
         private
 
         public initData_C

         interface initData_C
            module procedure initData_1d_C
            module procedure initData_2d_C
            module procedure initData_3d_C
            module procedure initData_4d_C
         end interface initData_C

         contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! initialization routines
!
            subroutine initData_1d_C(data,offt_1,offt_2)
               integer :: i, offt_1,offt_2
               complex(kind=kcomplex) :: data(:)

! aspg: for some reason, initializing the data here doesnt work starting
!       from offt_1 up to offt_2, but starting from 1 to the size
!       of the array works. wonder if this is a compiler issue or me
!       just making a mess..

               do i=1,(offt_2-offt_1)+1
                  data(i) = 0.0D0
               enddo

            end subroutine initData_1d_C


            subroutine initData_2d_C(data,offt1_i,offt2_i,offt1_j,offt2_j)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               complex(kind=kcomplex) :: data(:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     data(i,j) = 0.0D0
                  enddo
               enddo

            end subroutine initData_2d_C


            subroutine initData_3d_C(data,offt1_i,offt2_i,offt1_j,offt2_j,offt1_k,offt2_k)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               complex(kind=kcomplex) :: data(:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        data(i,j,k) = 0.0D0
                     enddo
                  enddo
               enddo

            end subroutine initData_3d_C


            subroutine initData_4d_C(data,offt1_i,offt2_i, &
     &                                      offt1_j,offt2_j, &
     &                                      offt1_k,offt2_k, &
     &                                      offt1_l,offt2_l )
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               integer :: l, offt1_l,offt2_l
               complex(kind=kcomplex) :: data(:,:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        do l=1,(offt2_l-offt1_l)+1
                           data(i,j,k,l) = 0.0D0
                        enddo
                     enddo
                  enddo
               enddo

            end subroutine initData_4d_C


      end module allocator_internal_init_C
!*
!*
!* Copyright (c) 2008-2013 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*
!
!

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module allocator_internal_init_R
      
         use allocator_parameters

         implicit none
         
         private
 
         public initData_R

         interface initData_R
            module procedure initData_1d_R
            module procedure initData_2d_R
            module procedure initData_3d_R
            module procedure initData_4d_R
         end interface initData_R

         contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! initialization routines
!
            subroutine initData_1d_R(data,offt_1,offt_2)
               integer :: i, offt_1,offt_2
               real(kind=kreal) :: data(:)

! aspg: for some reason, initializing the data here doesnt work starting
!       from offt_1 up to offt_2, but starting from 1 to the size
!       of the array works. wonder if this is a compiler issue or me
!       just making a mess..

               do i=1,(offt_2-offt_1)+1
                  data(i) = 0.0D0
               enddo

            end subroutine initData_1d_R


            subroutine initData_2d_R(data,offt1_i,offt2_i,offt1_j,offt2_j)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               real(kind=kreal) :: data(:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     data(i,j) = 0.0D0
                  enddo
               enddo

            end subroutine initData_2d_R


            subroutine initData_3d_R(data,offt1_i,offt2_i,offt1_j,offt2_j,offt1_k,offt2_k)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               real(kind=kreal) :: data(:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        data(i,j,k) = 0.0D0
                     enddo
                  enddo
               enddo

            end subroutine initData_3d_R


            subroutine initData_4d_R(data,offt1_i,offt2_i, &
     &                                      offt1_j,offt2_j, &
     &                                      offt1_k,offt2_k, &
     &                                      offt1_l,offt2_l )
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               integer :: l, offt1_l,offt2_l
               real(kind=kreal) :: data(:,:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        do l=1,(offt2_l-offt1_l)+1
                           data(i,j,k,l) = 0.0D0
                        enddo
                     enddo
                  enddo
               enddo

            end subroutine initData_4d_R


      end module allocator_internal_init_R
!*
!*
!* Copyright (c) 2008-2013 Andre Severo Pereira Gomes <andre.gomes@univ-lille1.fr>
!* All rights reserved.
!*
!* Redistribution and use in source and binary forms, with or without
!* modification, are permitted provided that the following conditions
!* are met:
!* 1. Redistributions of source code must retain the above copyright
!*    notice, this list of conditions and the following disclaimer.
!* 2. Redistributions in binary form must reproduce the above copyright
!*    notice, this list of conditions and the following disclaimer in the
!*    documentation and/or other materials provided with the distribution.
!*
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
!* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
!* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
!* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
!* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!* SUCH DAMAGE.
!*
!*
!
!

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module allocator_internal_init_I
      
         use allocator_parameters

         implicit none
         
         private
 
         public initData_I

         interface initData_I
            module procedure initData_1d_I
            module procedure initData_2d_I
            module procedure initData_3d_I
            module procedure initData_4d_I
         end interface initData_I

         contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! initialization routines
!
            subroutine initData_1d_I(data,offt_1,offt_2)
               integer :: i, offt_1,offt_2
               integer(kind=kinteger) :: data(:)

! aspg: for some reason, initializing the data here doesnt work starting
!       from offt_1 up to offt_2, but starting from 1 to the size
!       of the array works. wonder if this is a compiler issue or me
!       just making a mess..

               do i=1,(offt_2-offt_1)+1
                  data(i) = 0.0D0
               enddo

            end subroutine initData_1d_I


            subroutine initData_2d_I(data,offt1_i,offt2_i,offt1_j,offt2_j)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer(kind=kinteger) :: data(:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     data(i,j) = 0.0D0
                  enddo
               enddo

            end subroutine initData_2d_I


            subroutine initData_3d_I(data,offt1_i,offt2_i,offt1_j,offt2_j,offt1_k,offt2_k)
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               integer(kind=kinteger) :: data(:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        data(i,j,k) = 0.0D0
                     enddo
                  enddo
               enddo

            end subroutine initData_3d_I


            subroutine initData_4d_I(data,offt1_i,offt2_i, &
     &                                      offt1_j,offt2_j, &
     &                                      offt1_k,offt2_k, &
     &                                      offt1_l,offt2_l )
               integer :: i, offt1_i,offt2_i
               integer :: j, offt1_j,offt2_j
               integer :: k, offt1_k,offt2_k
               integer :: l, offt1_l,offt2_l
               integer(kind=kinteger) :: data(:,:,:,:)

               do i=1,(offt2_i-offt1_i)+1
                  do j=1,(offt2_j-offt1_j)+1
                     do k=1,(offt2_k-offt1_k)+1
                        do l=1,(offt2_l-offt1_l)+1
                           data(i,j,k,l) = 0.0D0
                        enddo
                     enddo
                  enddo
               enddo

            end subroutine initData_4d_I


      end module allocator_internal_init_I

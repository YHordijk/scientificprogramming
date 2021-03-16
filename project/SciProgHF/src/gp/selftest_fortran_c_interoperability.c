#include <stdio.h>

/* 
 This C function is called from Fortran as subroutine 
 to verify corectness of Fortran-to-C/C++ and back 
 data passing based upon ISO_C_BINDING.

 Written by Miro Ilias, Nobember 2015
*/

/* xlf90 linker requires underscore */
#if defined VAR_XLF
void selftest_fortran_c_interoperability_ 
#else 
void selftest_fortran_c_interoperability 
#endif
(int *arr_size, int *int_c_int, long *int_c_long,    
long int_c_long_array[], double *real_c_double,      
double real_c_double_array[], char string_c_char[], int *string_len, 
long *int_c_long_ptr,  double *real_c_double_ptr) 
{
/* simple manipulation: flip entering data through auxiliary variables and return them back  */
int i;  char c;
int  int_c_int_saved=*int_c_int;      
    *int_c_int  = int_c_int_saved;  
long  int_c_long_saved=*int_c_long;           
      *int_c_long  = int_c_long_saved; 
      int_c_long_saved=*int_c_long_ptr;
       *int_c_long_ptr=int_c_long_saved;
double  real_c_double_saved=*real_c_double;   
       *real_c_double = real_c_double_saved; 
       real_c_double_saved=*real_c_double_ptr; 
         *real_c_double_ptr=real_c_double_saved;

 for (i=0;i<*arr_size;i++) {
   real_c_double_saved=real_c_double_array[i];real_c_double_array[i]=real_c_double_saved;
   int_c_long_saved=int_c_long_array[i];int_c_long_array[i]=int_c_long_saved; }

 for (i=0;i<*string_len-1;i++) {
   c=string_c_char[i];string_c_char[i]=c; }
}


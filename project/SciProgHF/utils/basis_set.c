/*
 *     Copyright (c) 2018 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://www.diracprogram.org
 */

/*
 *  Program to read and write relativistic basis sets
 *
 *  Lucas Visscher
 *  Vrije Universiteit Amsterdam
 *
 *  May 2003
 *
 *  Last modifications: Miro Ilias, April 2006
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#define MAX_DOCUMENTATION 200   /* Maximum number of documentation lines                    */
#define MAX_REFERENCE 10        /* Maximum number of lines with literature references       */
#define MAX_COMMENTS 10         /* Maximum number of lines with atom-specific comments      */
#define MAX_LVALUE  8           /* Maximum value for angular momentum quantum number        */
#define MAX_BLOCK   50          /* Maximum number of function blocks of a particular type   */
#define MAX_EXPONENTS 200       /* Maximum number of different exponents                    */
#define MAX_COEFFICIENTS 1500   /* Maximum number of different coefficients                 */
#define MAX_ATOMS  120           /* Maximum number of atoms for a certain basis set type     */

#define DALTON_NUMBER_FORMAT "%5i%5i%5i\n"
#define DALTON_EXPONENT_FORMAT "%16.9le"
#define DALTON_COEFFICIENT_FORMAT "%12.9lf"
#define DYALL_NUMBER_FORMAT "%i %i %i %i %i"
#define DYALL_EXPONENT_FORMAT "%i %le"
#define DYALL_COEFFICIENT_FORMAT "%le"


/* Global variables */

const char l_name[MAX_LVALUE+1]={'S','P','D','F','G','H','I','J','K'};

int iprnt;

/* We may have different types of exponent relations  */

typedef enum
{
    non_family,     /* Individually l-optimized sets (e.g. non-relativistic) */
    dual_family,    /* Constrained l-optimized sets (d is subset s, p is subset f, etc.) */
    full_family,    /* All l-values are based on the same exponent series */
    j_optimized     /* l+ and l- sets have different exponents */
}   BASIS_TYPE;


/* We distinguish different types of contraction  */

typedef enum
{
    uncontracted,   /* No contraction is used */
    general,        /* All exponents are used in each contracted function (e.g. ANO-type sets) */
    segmented,      /* A set of exponents is only used in one contraction */
    mixed           /* Combination of the types given above */
}   CONTRACTION;

/* Text lines are used to write reference information */

typedef char LINE[132];


/* General information about the basis set is stored in struct BASIS */

typedef struct
{
    char name[20];                              /* Every basis has a name - it's output file name */
    char input_filename[20];                    /* Input file name - this file is to be read */
    BASIS_TYPE type ;                           /* Type of basis set */
    CONTRACTION contraction;                    /* Type of contraction */
    LINE documentation[MAX_DOCUMENTATION];      /* Documentation */
    LINE reference[MAX_REFERENCE];              /* Literature references */
    int number_of_atoms;                        /* Number of elements for which this basis set is defined */
}   BASIS;

/* The information for each atom is stored in the struct ATOOM
   Here we use arrays of dimension [MAX_BLOCK][3][MAX_LVALUE+1][2]
   where the second dimension distinguishes between l-, l+, and l(spinfree) respectively.
   The last dimension distinguishes between large and small (or pseudo-large in case of spinfree) components */

typedef struct
{
    LINE comments[MAX_COMMENTS];                                /* Commentary lines */
    int number;                                                 /* Atomic number */
    char symbol[4];                                             /* Element symbol */
    int l_max;                                                  /* Maximum l-value */
    int n_block[3][MAX_LVALUE+1][2];                            /* Number of basis function blocks for a particular l-value */
    int n_exponents[MAX_BLOCK][3][MAX_LVALUE+1][2] ;            /* Number of exponents for each basis function block */
    int n_orbitals[MAX_BLOCK][3][MAX_LVALUE+1][2] ;             /* Number of contracted functions in each basis function block */
    double exponents[MAX_EXPONENTS] ;                           /* Exponents */
    double coefficients[MAX_COEFFICIENTS] ;                     /* Contraction coefficients */
    double *exponent_start[MAX_BLOCK][3][MAX_LVALUE+1][2] ;     /* Pointer to start of exponent list for a particular function block */
    double *coefficient_start[MAX_BLOCK][3][MAX_LVALUE+1][2] ;  /* Pointer to start of coefficient list for a particular function block */
    double *exponent_end;                                       /* Pointer to end of exponent list (last exponent + 1) */
    double *coefficient_end;                                    /* Pointer to end of coefficient list (last coefficient + 1) */
}   ATOOM;

/* Function prototypes */

int read_dyall_stuff (BASIS *basis_set, ATOOM atoms[]);
int read_dyall_atom (FILE *file, ATOOM *atom);
int read_dyall_additional (FILE *file, ATOOM *atom, int l_max, int l_skip, int n_extra[]);
void read_block (FILE *file, int n_exp, int n_orbs, double *exponent, double *coefficient);
int move_to_section (FILE *file, const char text[], int skip_rewind);
int end_of_section (LINE line, const char text[]);
int write_dalton_style (BASIS basis_set, ATOOM atoms[]);
int write_element (FILE *file, ATOOM atom);
void write_block (FILE *file, int n_exp, int n_orbs, double *exp_start, double *coef_start);
void write_comment_lines (FILE *file, LINE comment[], int number_of_lines);
void combine_blocks (ATOOM *atom);

void process_main_args(int argc,char *argv[], char outna[], char inpna[]); /* process arguments  */
void print_help(void);

/* ============================================================================================================== */

//int main(void)
int main(int argc,char *argv[])
{
    int ierr = 0;
    BASIS basis_set;
    ATOOM atoms[MAX_ATOMS];

    /* Process the arguments  */
    process_main_args(argc,argv, basis_set.name, basis_set.input_filename);

    if (iprnt >= 0) {
    printf("\nmain: The input file to be read is: %s\n",basis_set.input_filename);
    printf("main: The output file to be produced is: %s\n\n",basis_set.name);
    }
 //   exit(10); 
/*

*/

    /*set_element_name (element_name);*/

    ierr = read_dyall_stuff (&basis_set, atoms);
    if (ierr == 1) {
      printf("Error in fopen(basis_set->input_filename) !!! ");
      exit(10);
    }

    write_dalton_style (basis_set, &atoms[0]);
    printf("Error code %i \n",ierr);

    return 0;
}
/* ============================================================================================================== */

void process_main_args(int argc,char *argv[], char outn[], char inpn[])
{
  int c;
 /* Flag set by `--verbose'. */
 static int verbose_flag;

  int is_out = 0;

 if (argc == 1)
 {
  printf("One or two arguments must be given ! Error exit !\n\n");
  print_help();
  exit(-1);
 }

 iprnt = 0;

 while (1)
 {
  static struct option long_options[] =
          {
               /* These options set a flag. */
               {"verbose", no_argument,       &verbose_flag, 1},
               {"brief",   no_argument,       &verbose_flag, 0},
               /* These options don't set a flag.
                  We distinguish them by their indices. */
               {"help",     no_argument,      0, 'h'},
               {"append",  no_argument,       0, 'b'},
               {"input",   required_argument, 0, 'i'},
               {"output",  required_argument, 0, 'o'},
               {"print",   required_argument, 0, 'p'},
               {0, 0, 0, 0}
             };
           /* getopt_long stores the option index here. */
           int option_index = 0;

           c = getopt_long (argc, argv, "hbi:o:p:",
                long_options, &option_index);

           /* Detect the end of the options. */
           if (c == -1)
             break;

           switch (c)
             {
             case 0:
               /* If this option set a flag, do nothing else now. */
               if (long_options[option_index].flag != 0)
                 break;
               printf ("option %s", long_options[option_index].name);
               if (optarg)
                 printf (" with arg %s", optarg);
               printf ("\n");
               break;

             case 'h':
               puts ("option -h\n");
               break;

             case 'b':
               puts ("option -b\n");
               break;

             case 'i':
          //     printf ("option -i with value `%s'\n", optarg);
               strcpy(inpn, optarg);
              // printf("The input file to be read is: %s\n\n",inpn);
           
               break;

             case 'o':
              // printf ("option -o with value `%s'\n", optarg);
               strcpy(outn, optarg);
              // printf("The output file to be produced is: %s\n\n",outn);
               is_out = 1;
               break;

             case 'p':
               printf ("option -p with value `%s'\n", optarg);
               // iprnt = (int) optarg;
               iprnt = atoi (optarg);
               printf("...entered print level number = %d \n",iprnt);
               break;

             case '?':
               /* getopt_long already printed an error message. */
               break;
            default:
              abort ();
             }
         }

       /* Instead of reporting `--verbose'
          and `--brief' as they are encountered,
          we report the final status resulting from them. */
       if (verbose_flag)
         puts ("verbose flag is set");

       /* Print any remaining command line arguments (not options). */
       if (optind < argc)
         {
           printf ("non-option ARGV-elements: ");
           while (optind < argc)
             printf ("%s ", argv[optind++]);
           putchar ('\n');
         }

    if (!is_out)
    {
     strcpy(outn,inpn);
     strcat(outn,".processed");
     if (iprnt > 4) {
     printf("The output file name is: %s\n\n",outn); }
    }
    return;
//  exit (0);
}

/* ============================================================================================================== */
   void print_help(void)
{
   printf("Usage: <this_program_exe> <parameter_list>)\n");
   printf(" ...where parameter list :                 \n");
   printf(" -i <basis_set_input_file>                 \n");
   printf(" -o <dirac_basis_set_output_file>          \n");
   printf(" -p <print_level (1-5) >                   \n");
//   printf(" -h                                        \n");

   return;
}
/* ============================================================================================================== */

int read_dyall_stuff (BASIS *basis_set, ATOOM atoms[])

/* Function to read the basis from the archive file provided by K.G. Dyall
   Return values :  0 - Normal completion
                    1 - Error while reading file
 */


{
    FILE *dyall_file;
    LINE line;
    int i, j, l, element;
    int l_max, n_extra[MAX_LVALUE+1];

    /* Open the file - the input file name is specified in the main routine */
   //  dyall_file = fopen("Dyall.qz", "r");
    dyall_file = fopen(basis_set->input_filename, "r");

    if (dyall_file == NULL) return 1;

    /* Initialize the fixed elements for this basis */
    // it's already initialized in the main routine...
    // strcpy(basis_set->name,"Dyall.VQZ");

    basis_set->type = non_family;
    basis_set->contraction = general;
    basis_set->number_of_atoms = 0;

    /* Read the reference lines */
    move_to_section (dyall_file,"REFERENCE",0);
    i = 0;
    do
    {
        fgets(line, 132, dyall_file);
        if (end_of_section (line,"REFERENCE") != 0) break;
        strcpy(basis_set->reference[i],line);
        i++;
    } while (i < MAX_REFERENCE);
    strcpy(basis_set->reference[i],"\0");

    /* Read the documentation lines */
    move_to_section (dyall_file,"DOCUMENTATION",0);
    i = 0;
    do
    {
        fgets(line, 132, dyall_file);
        if (end_of_section (line,"DOCUMENTATION") != 0) break;
        strcpy(basis_set->documentation[i],line);
        i++;
    } while ( i < MAX_DOCUMENTATION);
    strcpy(basis_set->documentation[i],"\0");

    /* Initialize the atoms */
    basis_set->number_of_atoms = 0;
    element = 0;

    /* Read the Hartree-Fock exponents and contraction coefficients */
    move_to_section (dyall_file,"RELATIVISTIC",0);
    do
    {
        fgets(line, 132, dyall_file);
        if (memcmp(line," **",3) == 0)
        {
            /* Write this also to output */
            printf ("Reading Hartree-Fock basis for %s",line);

            /* The first line is treated as the comment line for this element */
            strcpy (atoms[element].comments[0], line);
            strcpy (atoms[element].comments[1],"\0");

            /* Extract the element symbol for later use */
            sscanf (&line[3],"%s",&atoms[element].symbol);

            /* The first line should be the atomic number */
            fscanf(dyall_file,"%i",&atoms[element].number);

            read_dyall_atom (dyall_file, &atoms[element]);
            element++;
            basis_set->number_of_atoms++;
        };
    } while ( end_of_section (line,"RELATIVISTIC") == 0 && feof(dyall_file) == 0);

    /* Read the outer core correlation exponents */
    move_to_section (dyall_file,"OUTCOR_RELATIVISTIC",0);
    fscanf (dyall_file,"%i",&l_max);

    /* Read the number of extra functions for each l value */
    for  (l=0; l<=l_max; l++)
        fscanf (dyall_file,"%i",&n_extra[l]);

    /* Start loop over lines with additional exponents */
    fscanf (dyall_file,"%s",&line);
    while ( end_of_section (line,"OUTCOR_RELATIVISTIC") == 0 && feof(dyall_file) == 0)
    {
        /* Read the element symbol and find the corresponding atom, no match : store in new atom */
        printf ("Outer core correlation functions for %s \n",line);
        for (element=0; element<=basis_set->number_of_atoms; element++)
            if (strcmp(atoms[element].symbol,line) == 0) break;

        read_dyall_additional (dyall_file, &atoms[element], l_max, 0, n_extra);

        /* Prepare for the next element */
        fscanf (dyall_file,"%s",&line);
    }

    /* Read two blocks of valence correlating exponents */
    for (i=0; i<2; i++)
    {
    /* Read the valence correlating exponents */
    move_to_section (dyall_file,"VALCORR_RELATIVISTIC",i);
    fscanf (dyall_file,"%i",&l_max);

    /* Read the number of extra functions for each l value */
    for  (l=0; l<=l_max; l++)
        fscanf (dyall_file,"%i",&n_extra[l]);

    /* Start loop over lines with additional exponents */
    fscanf (dyall_file,"%s",&line);
    while ( end_of_section (line,"VALCORR_RELATIVISTIC") == 0 && feof(dyall_file) == 0)
    {
        /* Read the element symbol and find the corresponding atom, no match : store in new atom */
        printf ("Correlation functions for %s \n",line);
        for (element=0; element<=basis_set->number_of_atoms; element++)
            if (strcmp(atoms[element].symbol,line) == 0) break;

        read_dyall_additional (dyall_file, &atoms[element], l_max, 3, n_extra);

        /* Prepare for the next element */
        fscanf (dyall_file,"%s",&line);
    }
    }

    /* Read the valence polarization exponents */
    move_to_section (dyall_file,"VALPOL_RELATIVISTIC",0);
    fscanf (dyall_file,"%i",&l_max);

    /* Read the number of extra functions for each l value */
    for  (l=0; l<=l_max; l++)
        fscanf (dyall_file,"%i",&n_extra[l]);

    /* Start loop over lines with additional exponents */
    fscanf (dyall_file,"%s",&line);
    while ( end_of_section (line,"VALPOL_RELATIVISTIC") == 0 && feof(dyall_file) == 0)
    {
        /* Read the element symbol and find the corresponding atom, no match : store in new atom */
        printf ("Polarization functions for %s \n",line);
        for (element=0; element<=basis_set->number_of_atoms; element++)
            if (strcmp(atoms[element].symbol,line) == 0) break;

        read_dyall_additional (dyall_file, &atoms[element], l_max, 0, n_extra);

        /* Prepare for the next element */
        fscanf (dyall_file,"%s",&line);
    }

    return 0;
}
/* ============================================================================================================== */

int read_dyall_atom (FILE *file, ATOOM *atom)

/* Function to read one element on a basis set of Dyall archive (GRASP) type
   Return values :  0 - No errors
                    1 - Error reading basis set
 */

{
    int component, l, lj, n_exp, n_orbs, l_max;
    double *exp, *coef;
    int block;
    LINE line;

    /* Initialize the maximum l-value */
    l_max = -1;

    /* Pointer to start of exponents array */
    exp = &(atom->exponents[0]);
    /* Pointer to start of contraction coefficients array */
    coef = &(atom->coefficients[0]);

    /* Initialize block counter */
    for (component=0; component<2; component++)
        for (l=0; l<=MAX_LVALUE; l++)
            for (lj=0; lj<2; lj++)
                atom->n_block[lj][l][component] = 0;

    /* Loop over the blocks with exponents and coefficients */
    do
    {
        fgets(line,132,file);
        if (memcmp(line,"$BLOCK",6) == 0)
        {
            /* We found a block of exponents and coefficients,
               find out which one it is */

            fscanf(file,DYALL_NUMBER_FORMAT,
                   &component,&l,&lj,&n_exp,&n_orbs);
	    if (l > l_max) l_max = l;

            /* Store this information */
            block = atom->n_block[lj][l][component];
            atom->n_block[lj][l][component]++;
            atom->n_exponents[block][lj][l][component] = n_exp;
            atom->n_orbitals[block][lj][l][component] = n_orbs;
            atom->exponent_start[block][lj][l][component] = exp;
            atom->coefficient_start[block][lj][l][component] = coef;

            /* Read the exponents & coefficients */
            read_block(file, n_exp, n_orbs, exp, coef);

            /* Update the pointer */
            exp += n_exp;
            coef += n_orbs * n_exp;
        }
    }  while (memcmp(line,"$END_ELEMENT",11) != 0  && feof(file) == 0 );

    /* Keep information about maxima updated */
    atom->l_max = l_max;
    atom->exponent_end = exp;
    atom->coefficient_end = coef;

    return 0;
}
/* ============================================================================================================== */

void read_block (FILE *file, int n_exp, int n_orbs, double *exponent, double *coefficient)

/* Read block of functions that have the GRASP type format (function number, exponent, coefficients)  */

{
    int i, j, k;

    for (i=0; i<n_exp; i++)
    {

        fscanf(file,DYALL_EXPONENT_FORMAT, &k, exponent);   /* Read (a dummy and) the exponent */
        exponent++;                                         /* Point to the next exponent */

        for (j=0; j<n_orbs; j++)
        {
            fscanf(file, DYALL_COEFFICIENT_FORMAT, coefficient);    /* Read the coefficient */
            coefficient++;                                          /* Point to the next coefficient */
        }
    }
    return;
}
/* ============================================================================================================== */

int read_dyall_additional (FILE *file, ATOOM *atom, int l_max, int l_skip, int n_extra[])

/* Read additional (polarization/correlating functions) that are on one line */
/* l_skip indicates the l-value up to which functions should be skipped (they are already in the scf-section) */

{
    int i, j, l, lj, element, block;
    double additional[10], *exp, *coef;

    /* The l_max value will probably become larger */
    if (l_max > atom->l_max)
       atom->l_max = l_max;

    /* Read the exponents and make the contraction coefficients reflect the uncontracted situation */
    exp = atom->exponent_end;
    coef = atom->coefficient_end;

    for (l=0; l<=l_max; l++)
    {
        if (n_extra[l] == 0) continue;

        for (i=0; i<n_extra[l]; i++)
            fscanf(file,"%le",&additional[i]);

        /* Fix for uncontracted sets as combine_blocks will not detect identical exponents
           Do not take the polarization functions into account if they are already in the original set */
          if ( atom->n_block[0][l][0] > 0 &&  l <= l_skip ) continue;

        for (lj=0; lj<2; lj++)
        {
            if (l != 0 || lj == 1)
            {
                block = atom->n_block[lj][l][0];
                atom->n_exponents[block][lj][l][0] += n_extra[l];
                atom->n_orbitals[block][lj][l][0] += n_extra[l];

                atom->exponent_start[block][lj][l][0]= exp;
                atom->coefficient_start[block][lj][l][0]= coef;

                for (i=0; i<n_extra[l]; i++)
                {
                    *exp = additional[i];
                    exp++;

                    for (j=0; j<n_extra[l]; j++)
                    {
                        if (j == i)
                           *coef = 1.0;
                        else
                           *coef = 0.0;
                        coef++;
                    }
                }

                atom->n_block[lj][l][0]++;
            }
        }
    }


    /* Adjust pointers */
    atom->exponent_end = exp;
    atom->coefficient_end = coef;
}

/* ============================================================================================================== */

int move_to_section (FILE *file, const char text[], int skip_rewind)

/* Function to find begin of section on file
   Return values :  0 - Section was found
                    1 - Section was not found
 */


{
    int found=1;
    LINE line, target_line;

    /* Define the line that we are looking for */

    strcpy (target_line, "$BEGIN_");
    strcat (target_line, text);

    if (skip_rewind == 0) rewind (file);

    while ( feof(file) == 0)
    {
        fgets(line, 132, file);
        found = memcmp(line,target_line,strlen(target_line));
        if (found == 0) break;
    }

    return found;
}

/* ============================================================================================================== */

int end_of_section (LINE line, const char text[])

/* Function to check for end of section on file
   Return values :  0 - Line is not the end of the section
                    1 - End of section
 */

{
    LINE target_line;

    /* Define the line that we are looking for */

    strcpy (target_line, "$END_");
    strcat (target_line, text);

    if (memcmp(line,target_line,strlen(target_line)) == 0)
        return 1;
    else
        return 0;
}


/* ============================================================================================================== */

int write_dalton_style (BASIS basis_set, ATOOM atoms[])

/* Function to write DALTON-style library file
   Return values :  0 - Normal completion
                    1 - Error while writing file
                    2 - Basis set incompatible with DALTON-style library
 */


{
    FILE *basis_set_file;
    LINE supported[4];
    int element, i_line, n_lines, return_code, position;

    /* Open the file, the name of the set should be also the name of the file */

    basis_set_file = fopen(basis_set.name, "w");
    if (basis_set_file == NULL) return 1;

    /* Write the identifier */

    fprintf(basis_set_file, "$ %s\n", basis_set.name) ;

    /* Write the documentation and reference lines */
    fprintf(basis_set_file,"$\n$ Description\n$\n");
    write_comment_lines (basis_set_file, basis_set.documentation, MAX_DOCUMENTATION);
    fprintf(basis_set_file,"$\n$  REFERENCE\n");
    write_comment_lines (basis_set_file, basis_set.reference, MAX_REFERENCE);

    /* Loop over the active atoms */

    n_lines = 0;
    position = 0;
    for (element=0; element < basis_set.number_of_atoms; element++)
    {
        combine_blocks (&atoms[element]);                                           /* Make one block for each type */
        return_code = write_element (basis_set_file, atoms[element]);               /* Write the basis set for this atom */
        if (return_code != 0) return return_code;                                   /* Check for write or other errors */
        sprintf (&supported[n_lines][position]," %s",atoms[element].symbol);        /* Add to line(s) with supported elements (max 40 on each line) */
        position += strlen(atoms[element].symbol)+1;
        if (position > 128)
	{
	    n_lines++;
            position = 0;
        }
    }

    /* Write supported elements */

    fprintf(basis_set_file, "$Elements supported\n");
    for (i_line=0; i_line< n_lines + 1; i_line++)
        fprintf(basis_set_file, "$%s\n", supported[i_line]);

    return 0;
}

/* ============================================================================================================== */

void write_comment_lines (FILE *file, LINE comment[], int number_of_lines)
{
    int line = 0;

    while ( line < number_of_lines && strlen(comment[line]) != 0 )
    {
        fprintf(file, "$ %s", comment[line]) ;
        line++;
    }
}
/* ============================================================================================================== */

int write_element (FILE *file, ATOOM atom)

/* Function to write one atom, DALTON style
   Return values :  0 - Normal completion
                    1 - Error while writing file                          (Not implemented yet)/
                    2 - Basis set incompatible with DALTON-style output   (Not implemented yet)
 */


{
    const int large=0;  /* We only print large component functions */
    int block;        /* We always take only the first block (basis set should not be subdivided) */

    int l_value, ljval, n_exp, n_orbs;
    double *exp_start, *coef_start;

    /* We first need to identify this atom by its atomic number */
    fprintf(file, "a %i\n", atom.number) ;

    /* There may be some specific comments concerning this atom */
    write_comment_lines (file, atom.comments, MAX_COMMENTS);

    /* Loop over the l-values */

    for (l_value=0; l_value <= atom.l_max; l_value++)
    {
        fprintf(file, "$ %c-TYPE FUNCTIONS\n", l_name[l_value]) ;
        if (l_value > 0)
            ljval = 0;
        else
            ljval = 1;

        /* Note that we write only the l-1/2 blocks for l>0 ! */

        for (block=0; block<atom.n_block[ljval][l_value][large]; block++)
        {

            /* Find the size of these blocks */
            n_exp = atom.n_exponents[block][ljval][l_value][large];
            n_orbs = atom.n_orbitals[block][ljval][l_value][large];

            /* Initialize starting address */
            exp_start =  atom.exponent_start[block][ljval][l_value][large];
            coef_start =  atom.coefficient_start[block][ljval][l_value][large];

            /* Write this block to file */
            /* We do write the contraction coefficients but as they can not be used at present
               in DIRAC we write 0 instead of n_orbs */
            fprintf(file, DALTON_NUMBER_FORMAT, n_exp,0, 0);  /* Write the header */
            write_block(file, n_exp, n_orbs, exp_start, coef_start) ;

        }
    }

    return 0;
}

/* ============================================================================================================== */

void write_block (FILE *file, int n_exp, int n_orbs, double *exp_start, double *coef_start)

/* Writes block of exponents and coefficients in DIRAC/DALTON format */
/* Number of exponents should be written in calling function as different conventions are used */

{
    int i, j;
    double *exponent;       /* Initialize to the address of the first exponent of this block */
    double *coefficient;    /* Initialize to the address of the first coefficient of this block */

    exponent=exp_start;
    coefficient=coef_start;

    for (i=0; i<n_exp; i++)
    {
        fprintf(file, DALTON_EXPONENT_FORMAT, *exponent);   /* Write the exponent */
        exponent++;                                         /* Point to the next exponent */

        for (j=0; j<n_orbs; j++)
        {
            fprintf(file, DALTON_COEFFICIENT_FORMAT, *coefficient); /* Write the coefficient */
            coefficient++;                                          /* Point to the next coefficient */
            /* 6 coefficients should be written to line 1, 7 to the following lines */
            if ( (j+1)%7 == 0 && j+1<n_orbs ) fprintf(file, "\n");
        }
        fprintf(file, "\n");
    }
}


/* ============================================================================================================== */

void combine_blocks (ATOOM *atom)

/* Function to combine all functions of the same type into one block */
/* Does not recognize identical exponents yet */

{
    int component, l, lj, block, i, j;
    int nt_exp, nt_coef, n_exp, n_orbs, before, after;
    double *exp, *coef, exponents[MAX_EXPONENTS], coefficients[MAX_COEFFICIENTS];

    nt_exp = 0;
    nt_coef = 0;
    for (component=0; component<2; component++)
        for (l=0; l<= atom->l_max; l++)
            for (lj=0; lj<2; lj++)
                if (l > 0 || lj == 1)
                {
                    /* Count the total number of exponents of this type */
                    n_exp = 0;
                    n_orbs = 0;
                    for (block=0; block<atom->n_block[lj][l][component]; block++)
                    {
                        n_exp += atom->n_exponents[block][lj][l][component];
                        n_orbs += atom->n_orbitals[block][lj][l][component];
                    }

                    /* Copy the functions to the a temporary array */
                    before = 0;
                    after = n_orbs;
                    for (block=0; block<atom->n_block[lj][l][component]; block++)
                    {
                        exp = atom->exponent_start[block][lj][l][component];
                        coef = atom->coefficient_start[block][lj][l][component];
                        after -= atom->n_orbitals[block][lj][l][component];

                        for (i=0; i<atom->n_exponents[block][lj][l][component]; i++)
                        {
                            exponents[nt_exp]=*exp;
                            nt_exp++;
                            exp++;

                            /* For the coefficients we need to insert zeroes for coefficients that refer to exponents outside the block */
                            for (j=0; j<before; j++)
                            {
                                coefficients[nt_coef] = 0.0;
                                nt_coef++;
                            }

                            for (j=0; j<atom->n_orbitals[block][lj][l][component]; j++)
                            {
                                coefficients[nt_coef]=*coef;
                                nt_coef++;
                                coef++;
                            }

                            for (j=0; j<after; j++)
                            {
                                coefficients[nt_coef] = 0.0;
                                nt_coef++;
                            }

                        }

                        before += atom->n_orbitals[block][lj][l][component];
                     }

                    /* Overwrite the dimension parameters, not yet the exponents (we would overwrite blocks still to be copied !) */
                    atom->n_block[lj][l][component] = 1;
                    atom->n_exponents[0][lj][l][component] = n_exp;
                    atom->n_orbitals[0][lj][l][component] = n_orbs;
                 }


    /* Now we can overwrite the exponents and coefficients */
    exp = &(atom->exponents[0]);
    coef = &(atom->coefficients[0]);
    nt_exp = 0;
    nt_coef = 0;
    block = 0;

    for (component=0; component<2; component++)
        for (l=0; l<= atom->l_max; l++)
            for (lj=0; lj<2; lj++)
                if (l > 0 || lj == 1)
                {
                    n_exp = atom->n_exponents[block][lj][l][component];
                    n_orbs = atom->n_orbitals[block][lj][l][component];

                    /* Copy the functions from the temp array */
                    atom->exponent_start[block][lj][l][component] = exp;
                    atom->coefficient_start[block][lj][l][component] = coef;

                    for (i=0; i<n_exp; i++)
                    {
                        *exp = exponents[nt_exp];
                        nt_exp++;
                        exp++;

                        for (j=0; j<n_orbs; j++)
                        {
                            *coef = coefficients[nt_coef];
                            nt_coef++;
                            coef++;
                        }
                    }
                }

    /* Adjust pointers */
    atom->exponent_end = exp;
    atom->coefficient_end = coef;

    return;
}
/* ============================================================================================================== */


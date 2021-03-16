import sys


def contains_flag(sys_argv, flag):
    return (any(x for x in sys_argv if x.startswith('--{0}='.format(flag))))


def postprocess_args(sys_argv, arguments):

    # if --mpi is selected and compilers are not selected
    # then compilers default to mpif90, mpicc, and mpicxx
    if arguments['--mpi']:
        if not contains_flag(sys_argv, 'fc') and not contains_flag(sys_argv, 'cc') and not contains_flag(sys_argv, 'cxx'):
           arguments['--fc'] = 'mpif90'
           arguments['--cc'] = 'mpicc'
           arguments['--cxx'] = 'mpicxx'

    # if one of the compilers contains "mpi" and --mpi is not selected, it is probably a user error
    # in this case stop the configuration
    asking_for_mpi_compiler = False
    for flag in ['fc', 'cc', 'cxx']:
        if contains_flag(sys_argv, 'fc'):
            if 'mpi' in arguments['--fc']:
                asking_for_mpi_compiler = True
    if asking_for_mpi_compiler and not arguments['--mpi']:
        sys.stderr.write('ERROR: you ask for an MPI compiler but have not specified --mpi\n')
        sys.exit(1)

    return arguments

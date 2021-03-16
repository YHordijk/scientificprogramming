def configure(options, input_files, extra_args):
    """
    This function is used by runtest to configure runtest
    at runtime for DIRAC specific launch command and file naming.
    """

    from os import path
    from sys import platform

    launcher = 'pam'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    (inp, mol) = input_files

    if platform == "win32":
        exe = 'dirac.x.exe'
    else:
        exe = 'dirac.x'

    command = []
    command.append('python {0}'.format(launcher_full_path))
    command.append('--dirac={0}'.format(path.join(options.binary_dir, exe)))
    command.append('--noarch --nobackup')
    command.append('--inp={0} --mol={1}'.format(inp, mol))
    if extra_args is not None:
        command.append(extra_args)

    full_command = ' '.join(command)

    inp_no_suffix = path.splitext(inp)[0]
    mol_no_suffix = path.splitext(mol)[0]

    output_prefix = '{0}_{1}'.format(inp_no_suffix, mol_no_suffix)

    relative_reference_path = 'result'

    return launcher, full_command, output_prefix, relative_reference_path

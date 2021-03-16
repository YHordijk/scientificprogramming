
import os
import sys
import runtest_v1 as runtest


def write_stderr(log_file, s):
    """
    Writes s to stderr and to file log_file
    unless log_file is None.
    """
    if log_file:
        with open(log_file, 'w') as f:
            f.write(s)
    sys.stderr.write(s)


class Filter(runtest.Filter):

    def __init__(self):
        runtest.Filter.__init__(self)

    def add(self, *args, **kwargs):
        try:
            runtest.Filter.add(self, *args, **kwargs)
        except runtest.FilterKeywordError as e:
            sys.stderr.write(str(e))  # FIXME currently not written to any log file
            sys.exit(-1)


class TestRun(runtest.TestRun):

    def __init__(self, _file, argv):
        runtest.TestRun.__init__(self, _file, argv)
        self.return_code = 0

    def run(self, inp_files, mol_files, f=None, args='', accepted_errors=[], print_args=False):

        launch_script = os.path.normpath(os.path.join(self.binary_dir, 'pam'))
        if self.skip_run:
            sys.stdout.write('\nskipping actual run\n')
        else:
            if not os.path.exists(launch_script):
                sys.stderr.write('ERROR: launch script %s not found\n' % launch_script)
                sys.stderr.write('       have you set the correct --binary-dir (or -b)?\n')
                sys.stderr.write('       try also --help\n')
                sys.exit(-1)

        if sys.platform == "win32":
            dirac_exe = 'dirac.x.exe'
        else:
            dirac_exe = 'dirac.x'

        launcher = 'python "%s" --dirac=%s --noarch --nobackup %s' % (launch_script, os.path.join(self.binary_dir, dirac_exe), args)

        for inp in inp_files:
            inp_no_suffix = os.path.splitext(inp)[0]
            for mol in mol_files:
                mol_no_suffix = os.path.splitext(mol)[0]
                out = '%s_%s.out' % (inp_no_suffix, mol_no_suffix)
                if print_args:
                    sys.stdout.write('\nrunning test: %s %s with arguments=%s\n' % (inp_no_suffix, mol_no_suffix, args))
                else:
                    sys.stdout.write('\nrunning test: %s %s\n' % (inp_no_suffix, mol_no_suffix))
                
                command = launcher + ' --inp=%s --mol=%s' % (inp, mol)
                try:
                    runtest.TestRun.execute(self,
                                            command=command,
                                            accepted_errors=accepted_errors)
                    if f is None:
                        sys.stdout.write('finished (no reference)\n')
                    else:
                        try:
                            f.check(self.work_dir, '%s' % out, 'result/%s' % out, self.verbose)
                            sys.stdout.write('passed\n')
                        except IOError as e:
                            write_stderr(self.log, 'ERROR: could not open file %s\n' % e.filename)
                            sys.exit(-1)
                        except runtest.TestFailedError as e:
                            write_stderr(self.log, str(e))
                            self.return_code += 1
                        except runtest.BadFilterError as e:
                            write_stderr(self.log, str(e))
                            sys.exit(-1)
                        except runtest.FilterKeywordError as e:
                            write_stderr(self.log, str(e))
                            sys.exit(-1)
                except runtest.AcceptedError as e:
                    sys.stdout.write(str(e))
                except runtest.SubprocessError as e:
                    write_stderr(self.log, str(e))
                    sys.exit(-1)

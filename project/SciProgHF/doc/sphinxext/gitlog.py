# coding=utf-8

"""
This extension adds git info about page to context.
"""

import os
import re
import subprocess

from datetime import datetime

DEFAULT_DATE_FORMAT = '%d.%m.%Y %H:%M:%S'
LOG_FORMAT = 'name=%an email=%ae date=%at sha=%H'
LOG_REGEX = re.compile(u'^name=(?P<name>[\w\s]+) email=(?P<email>[\w@\.\s\-]+) date=(?P<date>\d{10}) sha=(?P<sha>\w{40})$', re.UNICODE)

_git_cache = {}


def git_installed():
    """
    Checks whether git binary is available.
    """
    if '_git_installed' not in _git_cache:
        process = subprocess.Popen('git --version', shell=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.communicate()
        _git_cache['_git_installed'] = (process.returncode == 0)
    return _git_cache['_git_installed']


def get_git_command(page_path, log_format):
    """
    Returns git log command string for given log format and page path.
    """
    return 'git log -1 --format="{0}" -- {1}'.format(log_format, page_path)


def get_git_root_path():
    """
    Returns root path of the git repository.
    """
    if '_git_root_path' not in _git_cache:
        root = check_subprocess_output('git rev-parse --show-toplevel', shell=True)
        _git_cache['_git_root_path'] = root.strip() if root else None
    return _git_cache['_git_root_path']


def get_page_path(cfg, root_path, pagename):
    """
    Returns absolute path to page by its pagename.
    """
    return os.path.join(root_path, cfg.gitlog_source_path, '{0}.rst'.format(pagename))


def get_git_log(cmd):
    """
    Runs the given command and returns its output.
    """
    return check_subprocess_output(cmd, shell=True)


def parse_log_info(info):
    """
    Parses the output of the git log command.
    """
    match = LOG_REGEX.search(info.strip())
    return match.groupdict() if match else None


def get_git_info(cfg, pagename):
    """
    Returns git info for given pagename.
    """
    root_path = get_git_root_path()
    page_path = get_page_path(cfg, root_path, pagename)
    cmd = get_git_command(page_path, LOG_FORMAT)
    info = get_git_log(cmd)
    return parse_log_info(info) if info else None


def format_date_string(date_str, date_format):
    """
    Returns formatted date string using given date format.
    """
    try:
        dt = datetime.fromtimestamp(float(date_str))
    except ValueError:
        return None
    return dt.strftime(date_format)


def fill_context(cfg, context, info):
    """
    Adds git info to page context.
    """
    context['git_author_name'] = info.get('name', '')
    context['git_author_email'] = info.get('email', '')
    context['git_commit_sha'] = info.get('sha', '')
    context['git_commit_date'] = format_date_string(info.get('date', ''), cfg.gitlog_date_format)


def page_context_handler(app, pagename, templatename, context, doctree):
    """
    Sets git info about page to page context.
    """
    if git_installed():
        info = get_git_info(app.config, pagename)

        if info:
            fill_context(app.config, context, info)


def setup(app):
    app.add_config_value('gitlog_source_path', DEFAULT_DATE_FORMAT, 'html')
    app.add_config_value('gitlog_date_format', DEFAULT_DATE_FORMAT, 'html')
    app.connect('html-page-context', page_context_handler)


def check_subprocess_output(*popenargs, **kwargs):
    """
    Runs command with arguments and return its output as a byte string. (Python 2.6, 2.7 compat)
    """
    if hasattr(subprocess, 'check_output'):
        return subprocess.check_output(*popenargs, **kwargs)
    return _check_output(*popenargs, **kwargs)


def _check_output(*popenargs, **kwargs):
    """
    Backport of Python 2.7 subprocess.check_output() for Python 2.6 support.
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

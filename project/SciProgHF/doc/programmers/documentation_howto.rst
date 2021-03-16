:orphan:
 

.. _documentation_howto:

How this documentation works and how you can contribute
=======================================================

This documentation is generated using `sphinx <http://sphinx.pocoo.org/>`_
which translates RST (reStructured Text) markup to html and latex.  Sphinx and
RST are extremely widely used in the python community.  The nice thing about
RST markup that it looks nice by itself in the terminal.  RST and Sphinx are
very simple and intuitive and for first steps consult the `reStructuredText
Primer <http://sphinx.pocoo.org/rest.html>`_ and the `Sphinx Markup Constructs
<http://sphinx.pocoo.org/markup/index.html>`_.  Typically RST markups is
auto-colored in the terminal so you can see when the markup is correct or
wrong.


Sphinx installation
-------------------

We recommend to install the required packages using virtualenv and not using
standard packages (you may need to install python-pip and python-virtualenv)::

  $ virtualenv venv

This will create a directory called ``venv`` in your current working directory.
You can create it in another place and you can give it a different name.
It may make sense to install it into ``~/venv``. It is no problem to have
20 different virtual environments on your computer.
Once you have set it up, you need to activate it::

  $ source venv/bin/activate

Next time you log on, you do not need to install it again, you just activate it
like above.  Once activated, you can install Python packages into the virtual
environment::

  $ pip install sphinx sphinx_rtd_theme sphinxcontrib-bibtex

The advantage of using virtualenv is that it is cross-platform,
you get the latest packages, you get an isolated and disposable
environment, and you do not require sudo permissions.

If you do not like the virtual environment anymore, you just delete it - it does
not pollute your system installation.

To learn more about virtual environments, see:
http://docs.python-guide.org/en/latest/dev/virtualenvs/.


Generation of html pages
------------------------

The central documentation web pages are generated automatically
from the RST files every hour.
But you can also generate the web pages by yourself. For this install
the sphinx-packages first (see above), then::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make -k html > build_html.log 2> build_html.log

Now point your browser to build/html/index.html. We recommend the developer to check the *build_html.log* file for possible errors
in generating the documentation.

How to add new sections and text
--------------------------------

Modify existing RST files or add new ones and reference them
in index.rst.

New theme
---------

The `sphinx_rtd_theme <https://github.com/snide/sphinx_rtd_theme>`_ has been customized and integrated as the new theme for Dirac. 
Theme extension can be installed using pip:

.. code-block:: bash

    pip install sphinx_rtd_theme

Here is example configuration with the fallback to the old theme:

.. code-block:: python

    try:
        import sphinx_rtd_theme
        html_theme = 'sphinx_rtd_theme'

        # set paths for loading local and also theme static resources
        html_theme_path = ['_themes', sphinx_rtd_theme.get_html_theme_path()]
        html_static_path = ['_static']

        # styles can be customized in this file
        html_context = {
            'css_files': ['_static/theme-custom.css'],
        }
    except ImportError:
        # fallback to old theme
        html_theme = 'dirac'
        html_theme_path = ['_themes']

If you want to customize the HTML templates, just copy the original ones from the theme and place them
inside the ``_templates`` folder inside the project. Sphinx will automatically prefer the ones in the ``_templates`` folder
over the global ones from the theme itself.

.. code-block:: bash

    # layouts can be customized in these paths
    _templates/layout.html
    _templates/footer.html

To override/customize the CSS styles, there is a css file placed under the ``_static`` folder called ``theme-custom.css``.

New format of citations
-----------------------

Citations have been updated to use the BibTeX format. This format is more maintainable and standardized. Support for this format is provided by the `sphinxcontrib-bibtex <http://sphinxcontrib-bibtex.readthedocs.org>`_ extension and it can be installed using pip:

.. code-block:: bash

    pip install sphinxcontrib-bibtex

Citations can be added using the new :rst:dir:`cite` directive:

.. code-block:: rest

    # old format
    [Dubillard2006]_

    # new format
    :cite:`Dubillard2006`

The new format also enables the option to download the citations in BibTeX.
This code will generate a link to static BibTeX file with all citations:

.. code-block:: rest

    :download:`BibTeX <zreferences.bib>`

The list of references can be added to the document using the
:rst:dir:`bibliography` directive:

.. code-block:: rest

    .. bibliography:: references.bib
      :enumtype: upperroman
      :style: diracstyle

.. note::

    Ensure that the bibliography directive is processed after all cites. The simplest way is to
    name the references page so that it starts with the "z" letter in order to be alphabeticaly last.
    The issue is being tracked `here <http://sphinxcontrib-bibtex.readthedocs.org/en/latest/usage.html#unresolved-citations-across-documents>`_.

When defining you own styles, please consider studying these `examples
<http://bazaar.launchpad.net/~pybtex-devs/pybtex/trunk/files/head:/pybtex/style/>`_,
as the documentation of this project is really poor. In order for citations to
be sorted alphabetically, the ``plain`` or ``alpha`` formatting must be used.

Available labels:

* ``alpha`` - Display data from citation fields.
* ``numbers`` - Displays numbers.

Available formatting:

* ``alpha`` - Uses ``alpha`` label, citations are sorted by ``author_year_title``.
* ``plain`` -  Uses ``number`` label, citations are sorted by ``author_year_title``.
* ``unsrt`` -  Uses ``number`` label, citations are sorted by order of appearance in documentation.
* ``unsrtalpha`` -  Uses ``alpha`` label, citations are sorted by order of appearance in documentation.

To customize citations style you can edit it in ``conf.py``, for example here we are using the citation key as the label:

.. code-block:: python

    from pybtex.plugin import register_plugin
    from pybtex.style.formatting.plain import Style
    from pybtex.style.labels.alpha import LabelStyle


    class DiracLabelStyle(LabelStyle):
        def format_label(self, entry):
            return entry.key


    class DiracStyle(Style):
        default_label_style = 'dirac'

    register_plugin('pybtex.style.labels', 'dirac', DiracLabelStyle)
    register_plugin('pybtex.style.formatting', 'diracstyle', DiracStyle)

Download mathjax formulas in LaTeX format:
------------------------------------------

The mathjax extension has been extended and renamed to **mathjaxplus**. The new **mathjaxplus** extension provides a way to download formulas in LaTeX format. There are also configuration parameters exposed to customize the generation process:

.. code-block:: python

    # disable or enable
    mathjax_generate_latex = True

    # style of save as link
    mathjax_latex_style = 'font-size:14px;text-align:right;'

    # title of the save as link
    mathjax_latex_title = 'Save as Latex'

Gitlog - display the GIT commit history for pages
-------------------------------------------------

A new extension called **gitlog** has been developed and integrated to the documentation. The motivation for this extension was to display information about the last git commit for each document in the documentation. There are some parameters exposed to configure the process:

.. code-block:: python

    # path to the git repository
    gitlog_source_path = '../doc/'

    # date format used to format git commit date in templates
    gitlog_date_format = '%d %b %Y, %H:%M:%S'

The extension adds following data available to the context of given document:

.. code-block:: bash

    git_author_name  # name of the commit author
    git_author_email # email address of the commit author
    git_commit_sha   # sha1 hash of the commit
    git_commit_date  # date of the commit

This code shows a simple usage in Jinja2 template:

.. code-block:: jinja

    {% if git_author_name %}
        Committer name: {{ git_author_name }}
        Committer email: {{ git_author_email }}
        Commit date: {{ git_commit_date }}
        Commit sha1 hash: {{ git_commit_sha }}
    {%- endif %}


Continuous integration and deployments
--------------------------------------

The current trend is to test your code on each push to the versioning system. We use this to ensure the stability
of the codebase and get valuable early feedback on our changes. There are some free CI services available which suits our needs.

Currently supported free CI services:

* `MagnumCI <https://magnum-ci.com>`_ 
* `Semaphore <https://magnum-ci.com>`_ 
* `Shippable <https://app.shippable.com>`_
* `Wercker <http://wercker.com/>`_
* `CircleCI <https://circleci.com/>`_
* `Codeship <https://codeship.com/>`_

Most of the CI services also offer ability to deploy built code. We use this feature to deploy the documentation after successful builds.
To deploy the code we use simple python script which wraps ``rsync`` command and adds some additional metadata to each deployed HTML. After
build, most of the services leaves the build environment untouched and so we can just deploy the built HTML files. Some services do not, and
we need to bootstrap the environment with system and python dependencies. To include build logs we need to rebuild documentation and capture
the standard and error outputs.

Example documentation build and deploy:

.. code-block:: bash

    BUILD_DIR=build_ci_name
    cd $BUILD_DIR
    SPHINX_LOG="sphinx-log.txt"
    DOXYGEN_LOG="doxygen-log.txt"
    SLIDES_LOG="slides-log.txt"
    make html > $SPHINX_LOG 2> $SPHINX_LOG
    make slides > $SLIDES_LOG 2> $SLIDES_LOG
    make doxygen > $DOXYGEN_LOG 2> $DOXYGEN_LOG
    cd ..
    python maintenance/deploy_doc.py --root=/path/to/documentation/root --user=user123 --host=my.hostname.com --port=22 --post_script=/path/to/post_deploy.py ci_name $BUILD_DIR

When deploying logs it is expected that logs are named according to this mapping:

.. code-block:: python

    mapping = {
        'sphinx': 'sphinx-log.txt',
        'slides': 'slides-log.txt',
        'doxygen': 'doxygen-log.txt',
    }

Deploy script is located at ``maintenance/deploy_doc.py`` and code is documented so feel free to study it when implementing new deployer.
For example this is implementation of :py:class:`Semaphore` deployer:

.. code-block:: python

    class Semaphore(Deployment):
        name = 'semaphoreci'

        def get_root_path(self):
            return os.environ['SEMAPHORE_PROJECT_DIR']

        def get_branch(self):
            return os.environ['BRANCH_NAME']

        def get_revision(self):
            return os.environ['REVISION']

As you can see we use environment variables to get information about the build. Each CI service provides its own set of environment variables so please refer to the respective CI documentation for more info. When implementing a new deployer you have to inherit from :py:class:`Deployment` and implement some methods to provide CI specific data about build:

.. code-block:: python

    def get_root_path(self):
        """
        Returns path to the root directory - git repository root.
        """
        raise NotImplementedError

    def get_revision(self):
        """
        Returns the current git revision.
        """
        raise NotImplementedError

    def get_branch(self):
        """
        Returns the current git branch.
        """
        raise NotImplementedError

After implementing the custom deployer class, you have to add it to CI mapping in order to be able to use it. You can do so in the :py:class:`main` function of the deployment script:

.. code-block:: python

    ci_mapping = {
        'semaphoreci': Semaphore,
        'codeship': Codeship,
        'magnumci': Magnum,
    }

To deliver the built HTML to the remote server we use ``rsync`` over SSH protocol. You can configure the script behaviour using keyword command line arguments:

--root         Remote root path, documentation will be synced to this path.
--user         Remote user.
--host         Remote hostname.
--port         Remote SSH port.
--post_script  Path to optional Python script on remote host.

After keyword arguments you have to provide two positional arguments.

- **ci** - CI name, as named in ``ci_mapping``.
- **build_dir** - The name of the build directory, for example ``build``, same as positional argument passed to the ``setup`` script.

Every operation on the remote server will be executed under the given user. If the ``post_script`` option is passed it will be executed after the rsync operation finishes. Example
post deploy script can be found in ``maintenance/post_deploy.py``. It used to generate nice overview of last documentation builds on the `main documentation site <http://dirac.umb.sk/doc/doc-deployments/>`_.

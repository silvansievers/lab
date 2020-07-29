Install Lab
-----------

Lab requires Python 3.6+ and Linux (e.g., Ubuntu). We recommend installing
Lab in a `Python virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_. This has the advantage
that there are no modifications to the system-wide configuration, and that
you can create multiple environments with different Lab versions (e.g.,
for different papers) without conflicts::

    # Install required packages, including virtualenv.
    sudo apt install python3 python3-venv

    # Create a new directory for your experiments.
    mkdir experiments-for-my-paper
    cd experiments-for-my-paper

    # If PYTHONPATH is set, unset it to obtain a clean environment.
    unset PYTHONPATH

    # Create and activate a Python 3 virtual environment for Lab.
    python3 -m venv --prompt my-paper .venv
    source .venv/bin/activate

    # Install Lab in the virtual environment.
    pip install -U pip wheel
    pip install lab  # or preferably a specific version with lab==x.y

    # Store installed packages and exact versions for reproducibility.
    # Ignore pkg-resources package (https://github.com/pypa/pip/issues/4022).
    pip freeze | grep -v "pkg-resources" > requirements.txt

If you want to install the latest development version and/or need to
change Lab itself, you can clone the Lab repo and install it in the
virtual environment::

    git clone https://github.com/aibasel/lab.git /path/to/lab
    pip install --editable /path/to/lab

The ``--editable`` flag installs the project in "editable mode", which
makes any changes under ``/path/to/lab`` appear immediately in the
installed package.

Please note that before running an experiment script you need to
activate the virtual environment with::

    source .venv/bin/activate

We recommend clearing the ``PYTHONPATH`` variable before activating the
virtual environment.

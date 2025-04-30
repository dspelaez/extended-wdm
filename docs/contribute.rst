Contributing guide
==================

We welcome contributions to the Extended Wavelet Directional Method (EWDM) project! Whether you're fixing bugs, adding new features, improving documentation, or sharing examples, your help is appreciated.

This guide provides instructions for setting up your development environment, running tests, and submitting contributions.

Setting up the development environment
--------------------------------------

1. **Clone the repository**

   First, clone the repository to your local machine:

   .. code-block:: bash

       git clone https://github.com/dspelaez/extended-wdm.git
       cd extended-wdm

2. **Create a virtual environment**

   It is recommended to use a virtual environment to manage dependencies:

   .. code-block:: bash

       conda create --name ewdm-dev python=3.8 # or higher
       conda activate ewdm-dev


3. **Install dependencies**

   Install the package in editable mode along with development and testing dependencies:

   .. code-block:: bash

       pip install -e '.[test]'

   This will install EWDM along with packages like `pytest` for testing.

Running tests
-------------

To ensure that your changes don't break existing functionality, run the test suite:

.. code-block:: bash

    pytest tests

This command will execute all tests located in the `tests/` directory.

Contributing new data sources
-----------------------------

The `ewdm.sources` module contains implementations for various data sources. If you have access to a new type of wave data (e.g., from a specific buoy or sensor), consider adding support for it:

1. **Implement a new source class**

   Create a new class in `ewdm/sources.py` that handles data loading and preprocessing for your source.

2. **Write tests**

   Add tests for your new source in the `tests/` directory to ensure its functionality.

3. **Update documentation**

   Document your new source in the appropriate section of the documentation to help others understand how to use it.

Adding examples to the gallery
------------------------------

We encourage contributors to add examples demonstrating the use of EWDM. To integrate your example into the documentation using Sphinx-Gallery, follow these steps:

1. **Create a Python script**

   Develop a Python script that showcases a specific use case or analysis using EWDM. Ensure that your script includes:

   - A docstring at the beginning, formatted in reStructuredText (reST), which will serve as the description in the gallery.

   - Code that generates plots or outputs, as these will be captured and displayed in the gallery.

   For guidance on structuring your script, refer to the `Sphinx-Gallery syntax documentation <https://sphinx-gallery.github.io/stable/syntax.html>`_.

2. **Place the script in the `docs/gallery/` directory**

   Save your script in the `docs/gallery/` directory of the repository. This directory is configured to be processed by Sphinx-Gallery during the documentation build.

3. **Build the documentation**

   After adding your script, build the documentation to ensure that your example appears correctly in the gallery:

   .. code-block:: bash

       make clean
       make html

   Navigate to the generated HTML files to verify that your example is displayed as intended.

By following these steps, your example will be automatically processed by Sphinx-Gallery and included in the documentation's gallery section. This allows users to view the example, download the script, and see the generated outputs directly within the documentation.


Adding standalone Jupyter notebooks
-----------------------------------

It is also possible to add standalone Jupyter Notebooks.

1. **Create a Jupyter Notebook**

   Develop a notebook that showcases a specific use case or analysis using EWDM.

2. **Add them to the `notebooks/` directory**

   Place your notebook in the `notebooks/` directory of the repository.

3. **Run them in the cloud**

    - You can run the Jupyter Notebooks in the cloud using `Binder <https://mybinder.org/v2/gh/dspelaez/extended-wdm/HEAD?labpath=notebooks>`_.
    - Wait for Binder to launch and then navigate to the desired `.ipynb` file in the Jupyter interface or create a new one.
    - Binder sessions are temporary, so download your notebooks to keep changes.


Contribution workflow
---------------------

To contribute to EWDM:

1. **Open an issue**

   Before making significant changes, open an issue to discuss your ideas with the maintainers.

2. **Fork the repository**

   Create a personal fork of the repository on GitHub.

3. **Create a feature branch**

   Develop your changes in a new branch:

   .. code-block:: bash

       git checkout -b feature/your-feature-name

4. **Testing**

   Write tests and make sure all tests are passed

   .. code-block:: bash

       pytest tests

5. **Commit and push changes**

   Commit your changes with clear messages and push to your fork:

   .. code-block:: bash

       git push origin feature/your-feature-name

6. **Submit a Pull Request**

   Open a pull request from your feature branch to the main repository's `main` branch.

7. **Review process**

   Engage in the review process by addressing feedback and making necessary adjustments.

Thank you for considering contributing to EWDM! Your efforts help improve the toolkit for the entire oceanographic community.

For more information, visit our documentation: https://extended-wdm.readthedocs.io/en/latest/

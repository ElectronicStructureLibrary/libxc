import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='pylibxc',
        version="4.0.1",
        description=
        'PyLibxc is a python-bound C library of exchange-correlation functionals for density-functional theory.',
        author='LibXC Authors',
        author_email='marques@tddft.org',
        url="https://gitlab.com/libxc/libxc",
        license='MPL',
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',    # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
        long_description="""
Libxc is a library of exchange-correlation functionals for
density-functional theory. The aim is to provide a portable, well
tested and reliable set of exchange and correlation functionals that
can be used by a variety of programs.

For more information, please check the manual at
http://www.tddft.org/programs/Libxc
""")

try:
    from setuptools import setup
except ImportError:
    try:
        import ez_setup
        ez_setup.use_setuptools()
    except ImportError:
        raise ImportError("wellFARE could not be installed, probably because"
            " neither setuptools nor ez_setup are installed on this computer."
            "\nInstall ez_setup ([sudo] pip install ez_setup) and try again.")

from setuptools import setup, find_packages


from setuptools import setup, find_packages

setup(name='wellfare',
      version='0.1.0',
      author='Valentin',
    description='',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="",
    scripts=['bin/wellfare'],
    install_requires= ['docopt', 'xlrd', 'numpy', 'scipy'],
    packages= find_packages(exclude=['docs','examples']))

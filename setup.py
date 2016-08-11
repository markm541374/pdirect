from setuptools import setup
from setuptools.extension import Extension



with open('pdirect/VERSION') as version_file:
    version = version_file.read().strip()

#def readme():
#    with open('pdirect/README.md') as f:
#        return f.read()

compile_flags = ['-O3']

from numpy import get_include as npgi
extensions = [
    Extension(name ="pdirect/core/direct",
              sources = ["pdirect/core/direct.c"],
              include_dirs = ['.','core',npgi()],
              extra_compile_args=compile_flags
    )
]

setup(name='pdirect',
      version=version,
      description='direct',

      url='https://github.com/markm541374/pdirect',
      author='markm541374',
      license='None',
      packages=['pdirect'],
      package_dir={'pdirect':'pdirect'},
      package_data={'pdirect':['VERSION','README.md']},
      ext_modules= extensions,
      zip_safe=False)

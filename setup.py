from setuptools import setup

setup(
    name='DoZen',
    version='0.1',
    description='Load, process, plot, and QC electromagnetic data from Zonge z3d files',
    url='https://git.multiphysics-mva.org/dozen.git',
    author='Andy McAliley',
    author_email='andy.mcaliley@gmail.com',
    packages=['dozen',],
    license='MIT license',
    # install_requires = ['pip','numpy','pandas','scipy','matplotlib','jupyter','pyviz==0.10','pyproj==1.9.6'],
    long_description=open('README.md').read(),
    zip_safe=False
)

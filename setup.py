from setuptools import setup

install_requires = ['pip','numpy','pandas','scipy','pyproj==1.9.6']
extras_requires = ['matplotlib','jupyter','pyviz==0.10']

setup(
    name='DoZen',
    version='0.1',
    description='Load, process, plot, and QC electromagnetic data from Zonge z3d files',
    url='https://github.com/AndyMcAliley/DoZen.git',
    author='Andy McAliley',
    author_email='andy.mcaliley@gmail.com',
    packages=['dozen',],
    license='MIT license',
    install_requires = install_requires,
    extras_requires = extras_requires,
    long_description=open('README.md').read(),
    zip_safe=False
)

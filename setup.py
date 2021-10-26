import re
from os import path
from setuptools import setup, find_packages


here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), 'r', encoding='utf-8') as fh:
    long_description = fh.read()


with open(path.join(here, 'protmapper', '__init__.py'), 'r') as fh:
    for line in fh.readlines():
        match = re.match(r'__version__ = \'(.+)\'', line)
        if match:
            version = match.groups()[0]
            break
    else:
        raise ValueError('Could not get version from protmapper/__init__.py')


def main():
    install_list = ['requests', 'boto3', 'pystow>=0.1.0']

    setup(name='protmapper',
          version=version,
          description='Map protein sites to human reference sequence.',
          long_description=long_description,
          long_description_content_type='text/markdown',
          author='John A. Bachman',
          author_email='john_bachman@hms.harvard.edu',
          url='https://github.com/indralab/protmapper',
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
          keywords=['protein', 'proteomics', 'sequence', 'alignment',
                    'assembly', 'post-translational', 'modification'],
          packages=find_packages(),
          install_requires=install_list,
          tests_require=['nose'],
          extras_require={'rest_api': ['flask', 'flask_cors']},
          include_package_data=True,
          entry_points={'console_scripts':
                        ['protmapper = protmapper.cli:main']},
        )


if __name__ == '__main__':
    main()

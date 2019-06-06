import sys
from setuptools import setup


def main():
    install_list = ['requests', 'rdflib', 'boto3']

    setup(name='protmapper',
          version='0.0.11',
          description='Map protein sites to human reference sequence.',
          long_description=('The protmapper is a tool to map inconsistent '
                            'protein sites (i.e., not matching the human '
                            'reference sequence) found in PTM databases and '
                            'the scientific literature to corresponding '
                            'positions on the human reference sequence.'),
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
          packages=['protmapper', 'protmapper.rest_api'],
          install_requires=install_list,
          tests_require=['nose'],
          extras_require={'rest_api': ['flask', 'flask_cors']},
          include_package_data=True,
        )


if __name__ == '__main__':
    main()

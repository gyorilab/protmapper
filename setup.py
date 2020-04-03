import os
from setuptools import setup


readme_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                           'README.md')
with open(readme_path, 'r', encoding='utf-8') as fh:
    long_description = fh.read()


def main():
    install_list = ['requests', 'rdflib', 'boto3']

    setup(name='protmapper',
          version='0.0.15',
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
          packages=['protmapper', 'protmapper.rest_api'],
          install_requires=install_list,
          tests_require=['nose'],
          extras_require={'rest_api': ['flask', 'flask_cors']},
          include_package_data=True,
          entry_points={'console_scripts':
                        ['protmapper = protmapper.cli:main']},
        )


if __name__ == '__main__':
    main()

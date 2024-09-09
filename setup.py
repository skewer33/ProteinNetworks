from setuptools import setup, find_packages

def readme():
  with open('README.md', 'r', encoding="utf-8") as f:
    return f.read()

setup(
  name='ProteinNetworks',
  version='0.1.2',
  author='Mokin Yakov',
  author_email='mokinyakov@mail.ru',
  description='Module for working with protein networks (gene ontology, enrichment, protein-protein interactions, etc.)',
  long_description=readme(),
  long_description_content_type='text/markdown',
  url='https://github.com/skewer33/ProteinNetworks.git',
  packages=find_packages(),
  install_requires=[
    'pandas>=2.0.1',
    'stringdb==0.1.5',
    'tabulate==0.9.0'],
  classifiers=[
    'Programming Language :: Python :: 3.12',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent'
  ],
  keywords='proteins interactions PPI networks enrichment STRINGdb Bioilogical-Processes Molecular-Functions Cellular-Components Gene-Ontology',
  project_urls={
    'Documentation': 'https://github.com/skewer33/ProteinNetworks/blob/main/README.md'
  },
  python_requires='>=3.7'
)

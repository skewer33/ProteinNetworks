from setuptools import setup, find_packages
import re
from pathlib import Path

def readme():
  with open('README.md', 'r', encoding="utf-8") as f:
    return f.read()

def get_version():
    init_path = Path(__file__).parent / 'ProteinNetworks' / '__init__.py'
    content = init_path.read_text()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Cannot find version in __init__.py")
  

def get_requirements():
    req_path = Path(__file__).parent / 'requirements.txt'
    with req_path.open() as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]
  
# prepare requirements and split heavy/compiled deps into extras to avoid mandatory builds
_reqs = get_requirements()
_heavy = {'matplotlib', 'leidenalg', 'igraph', 'umap-learn', 'python-igraph'}
_install_requires = [r for r in _reqs if r.split('==')[0] not in _heavy]
_extras = {'full': [r for r in _reqs if r.split('==')[0] in _heavy]}

setup(
  name='ProteinNetworks',
  version=get_version(),
  author='Mokin Yakov',
  author_email='mokinyakov@mail.ru',
  description='Module for working with protein networks (gene ontology, enrichment, protein-protein interactions, etc.)',
  long_description=readme(),
  long_description_content_type='text/markdown',
  url='https://github.com/skewer33/ProteinNetworks.git',
  packages=find_packages(),
  install_requires=_install_requires,
  extras_require=_extras,
  classifiers=[
    'Programming Language :: Python :: 3.12',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent'
  ],
  keywords='proteins interactions PPI networks enrichment STRINGdb Bioilogical-Processes Molecular-Functions Cellular-Components Gene-Ontology',
  project_urls={
    'Documentation': 'https://github.com/skewer33/ProteinNetworks/blob/main/README.md'
  },
  python_requires='>=3.7',
  include_package_data=True
)

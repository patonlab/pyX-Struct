from setuptools import setup
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'pyX-Struct',
  packages = ['pyX-Struct'],
  version = '1.0',
  description = 'Scrape Geometric X-ray Data from the Cambridge Structural Database ',
  long_description=long_description,
  long_description_content_type='text/markdown',
  author = 'Paton Research Group',
  author_email = 'robert.paton@colostate.edu',
  url = 'https://github.com/bobbypaton/pyX-Struct',
  download_url = 'https://github.com/bobbypaton/pyX-Struct/archive/1.0.zip',
  keywords = ['x-ray structure', 'CCDC', 'SMILES', 'python'],
  classifiers = [],
  install_requires=["numpy","seaborn","pandas","matplotlib"],
  python_requires='>=2.6',
  include_package_data=True,
)


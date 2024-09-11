from setuptools import setup, find_packages

setup(
    name='FastHareComposite',
    version='0.3.5',
    author='Thang N. Dinh',
    author_email='thangdn@gmail.com',
    description='A composite for using FastHare reduction with the samplers from D-Wave Ocean SDK',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/thang-dinh/FastHare",
    packages=find_packages(),
    install_requires=[
        'dimod',
        'fasthare',  
        'networkx'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)

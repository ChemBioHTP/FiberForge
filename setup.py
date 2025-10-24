from setuptools import setup, find_packages

setup(
    name='FiberForge',
    version='0.1.0',
    packages=find_packages(),  # Automatically discover and include all packages in the project

    # Project metadata
    author='Kieran Nehil-Puleo',
    author_email='nehilkieran@gmail.com',
    description='A package for building and characterizing amyloid fibrils using MD',
    url='https://github.com/kierannp/FiberForge',
    license='MIT',

    # Additional classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],

    # Other configurations
    zip_safe=False,
)

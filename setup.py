from setuptools import setup

setup(
    name = 'acheron',
    version = '0.1.0',
    packages = ['acheron'],
    entry_points = {
        'console_scripts': [
            'acheron = acheron.__main__:main'
        ]
    })

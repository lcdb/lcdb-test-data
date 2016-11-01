from setuptools import setup
setup(
    name='lcdb_test_data',
    author='Ryan Dale',
    description='Generates test data and provides access to it',
    license='MIT',
    include_package_data=True,
    packages=['lcdb_test_data'],
    entry_points={
        'console_scripts':
        [
            'build_example_data = lcdb_test_data.build:main',
        ],
    },
)

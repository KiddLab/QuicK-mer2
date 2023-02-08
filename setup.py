from setuptools import setup

setup(name='qm2_human_rarity_prj',
      version='0.1.0',
      packages=['qm2_human_rarity'],
      entry_points={
            'console_scripts': [
                  'test_run = qm2_human_rarity.test:main'
                  'test_dups = qm2_human_rarity.compare_against_1000:main'
            ]
      })

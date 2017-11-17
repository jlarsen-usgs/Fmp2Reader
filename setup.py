from distutils.core import setup

setup(name="fmp2reader",
      version="0.1",
      author="Joshua D Larsen",
      author_email="jlarsen@usgs.gov",
      packages=["fmp2reader"],
      install_requires=["pyshp", "numpy>=1.7"],
      description="""Python utility to read Fmp files from OWHM2 and create
                     shapefiles from that data, incomplete, contribute to development
                     through github!"""
      )

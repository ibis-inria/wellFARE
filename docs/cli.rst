Command-line interface
=======================

This is a copy of the help for the command line interface, which you can also get by typing ``wellfare -h`` in a terminal

::
    
    Usage:
      wellfare growth <volumefile> <outputfile>
      wellfare activity <volumefile> <fluofile> <dR> <outputfile>
      wellfare protein <volumefile> <fluofile> <dR> <dP> <outputfile>
      wellfare synchronize <reference> <smin=100> <smax=-100> <file1> <file2>...
      wellfare outliers



    Documentation:
    ===============

    WellFARE_cli is a command-line interface to the python package
    wellFARE, for use with other languages such as MATLAB, PHP, etc. 


    Currently supported are functions for the computation of growth 
    rate, promoter activity, associated protein concentration,
    outliers removal, and curves synchronization.

    wellfare-cli communicates with other languages by the intermediary
    of files. All data files read or written by wellfare will be of
    the following form

        t (min)   volume
        0.000     1.001
        3.000     1.100
        etc.

    (the first line is optional, separators are spaces)

    API
    ====

    growth

      Estimates the growth rate from volume data

      Parameters
      -----------

      volumefile

        Name of the file containing the volume time-series data.

      outputfile
        
        Name of the file to write the result in.



    promact
      
      Estimates (proportionally) the promoter activity (or
      synthesis rate) responsible for the observed fluorescence
      signal.
      Supposes that the production of GFP occurs in one step,
      and produces a proptional estimate of 

      Parameters
      -----------

        volumefile

          Name of the file containing the volume time-series data.

        fluofile

          Name of the file containing the fluorescence time-series data.

        dR

          Degradation rate of the reporter, in the same unit as the time
          in the data files.

        outputfile
        
          Name of the file to write the result in.



    protein

      Estimates (proportionally) the promoter activity (or
      synthesis rate) responsible for the observed fluorescence
      signal.
      Supposes that the production of GFP occurs in one step,
      and produces a proptional estimate of 

      Parameters
      -----------

        volumefile

          Name of the file containing the volume time-series data.

        fluofile

          Name of the file containing the fluorescence time-series data.

        dR

          Degradation rate of the reporter, in the same unit as the time
          in the data files.

        outputfile
        
          Name of the file to write the result in.
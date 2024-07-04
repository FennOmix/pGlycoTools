# Python scripts in pGlyco


### Using XIC.py to extract 15N/13C-labeled and 13C intensities

1. Change `glyco.ini` to `glyco-15N.ini` (`glycoini=glyco-15N.ini`) in pGlyco.cfg (for instance `C:/glyco/pGlyco.cfg`) in the pGlyco3/pGlycoNovo result folder.
2. Run `python XIC.py -p C:/glyco/pGlyco.cfg`, then the `PeakArea` column will contain the peak areas of 15N-labeled glycopeptides.

The step is similar for `glyco-13C.ini`.

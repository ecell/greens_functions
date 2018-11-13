all: 2DRadAbs R31P00

R31P00:
	platex GreensFunctionR31P00.tex && platex GreensFunctionR31P00.tex && dvipdfmx GreensFunctionR31P00.dvi

2DRadAbs:
	platex GreensFunction2DRadAbs.tex && platex GreensFunction2DRadAbs.tex && dvipdfmx GreensFunction2DRadAbs.dvi

.PHONY: clean
clean:
	rm *.aux *.dvi *.log *.pdf

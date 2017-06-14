LAT = pdflatex 
LATFLAGS = -shell-escape

0:
	 $(LAT) $(LATFLAGS) 0.tex

1:
	 $(LAT) $(LATFLAGS) 1.tex

clean:
	 rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz

all: 0 1 
LAT = pdflatex 
LATFLAGS = -shell-escape

1:
	 $(LAT) $(LATFLAGS) 1.tex

2:
	 $(LAT) $(LATFLAGS) 2.tex

clean:
	 rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz

all: 1 2

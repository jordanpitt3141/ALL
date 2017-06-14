LAT = pdflatex 
LATFLAGS = -shell-escape

0:
	 $(LAT) $(LATFLAGS) 0.tex

1:
	 $(LAT) $(LATFLAGS) 1.tex

2:
	 $(LAT) $(LATFLAGS) 2.tex

3:
	 $(LAT) $(LATFLAGS) 3.tex

4:
	 $(LAT) $(LATFLAGS) 4.tex

clean:
	 rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz

all: 0 1 2 3 4 
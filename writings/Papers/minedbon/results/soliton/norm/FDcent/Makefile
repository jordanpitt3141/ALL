LAT = pdflatex
LATFLAGS = -shell-escape

o2:
	$(LAT) $(LATFLAGS) sto2.tex
o2E:
	$(LAT) $(LATFLAGS) sto2E.tex

clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5

all: o2 


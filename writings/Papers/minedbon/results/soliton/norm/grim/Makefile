LAT = pdflatex
LATFLAGS = -shell-escape

o1:
	$(LAT) $(LATFLAGS) sto1.tex
o1E:
	$(LAT) $(LATFLAGS) sto1E.tex

clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5

all: o1


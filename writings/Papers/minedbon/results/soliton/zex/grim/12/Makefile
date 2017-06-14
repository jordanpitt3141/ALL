LAT = pdflatex
LATFLAGS = -shell-escape

u:
	$(LAT) $(LATFLAGS) plotsolu.tex
h:
	$(LAT) $(LATFLAGS) plotsolh.tex

clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5

all: u h


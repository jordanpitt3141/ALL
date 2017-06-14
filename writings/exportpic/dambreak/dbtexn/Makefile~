LAT = pdflatex
LATFLAGS = -shell-escape

both:
	$(LAT) $(LATFLAGS) both.tex

clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz

all: both


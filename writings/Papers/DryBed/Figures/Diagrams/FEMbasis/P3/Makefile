SRC_FILES := $(wildcard *.tex)
OUT_FILES :=  $(addsuffix .pdf,$(basename $(wildcard *.tex)))

$(OUT_FILES):
	pdflatex -shell-escape  $(addsuffix .tex,$(basename $@))

all: $(OUT_FILES)

clean:
	rm -f *.aux *.auxlock *.log *.dep *.dpth *.md5 *temp* *.pdf *.synctex.gz *.gnuplot *.table

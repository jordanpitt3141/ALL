SRC_DIR := .
SRC_FILES := $(wildcard $(SRC_DIR)/*.tex)

all: 
	pdflatex -shell-escape $(SRC_FILES)

clean:
	rm -f *.aux *.auxlock *.log *.dep *.dpth *.md5 *temp* *.pdf

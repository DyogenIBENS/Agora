
MERMAID_WEB_SERVICE = "https://mermaid.ink/img/"
MMD_FILES = $(wildcard *.mmd)
JPG_FILES = $(MMD_FILES:.mmd=.jpg)

all: HowTo.pdf
howto: HowTo.pdf
diagrams: $(JPG_FILES)

HowTo.pdf: HowTo.md diagrams
	pandoc HowTo.md -o HowTo.docx
	lowriter --convert-to pdf HowTo.docx

%.jpg: %.mmd
	$(eval BASE64=$(shell base64 < "$?" | tr -d '\n'))
	curl "$(MERMAID_WEB_SERVICE)$(BASE64)" -o $@

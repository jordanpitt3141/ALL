from subprocess import call

namelist = ["*.tex","*.tikz","*.log","*.aux","*.auxlock","*.dpth","*.md5"]

for extension in namelist:
    call(["find",".","-name",extension,"-type","f","-delete"])

#find .-name   "*.auxlock" -type f -delete

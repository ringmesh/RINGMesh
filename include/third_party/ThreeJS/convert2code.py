import os
import glob
from os.path import basename

prefix = "out << \""
suffix = "\" << std::endl ;"

Fichier = open("output_code.txt",'w')
files = glob.glob("*.js")
for file in files:
	fileName, fileExtension = os.path.splitext(file)
	fileName = fileName.replace('.', '_')
	newFile = fileName + ".txt"
	Fichier.write( "void print_" + fileName + "( std::ofstream& out ) { \n" )
	with open(file, "r") as ins:
		for line in ins:
			line = line.rstrip('\n')
			line = line.replace('\\', '\\\\')
			line = line.replace('"', '\\"')
			Fichier.write( '%s%s%s\n' % (prefix, line, suffix) )
	Fichier.write( "}\n\n" )
Fichier.close()


Fichier = open("main.txt",'w')
Fichier.write( "void print_main( std::ofstream& out ) { \n" )
with open("main.html", "r") as ins:
	for line in ins:
		line = line.rstrip('\n')
		line = line.replace('\\', '\\\\')
		line = line.replace('"', '\\"')
		Fichier.write( '%s%s%s\n' % (prefix, line, suffix) )
Fichier.write( "}\n\n" )
Fichier.close()


OBJS = objs/Vector3.o objs/Matrix.o objs/Quaternion.o objs/BandContrast.o objs/BandContrastWriter.o objs/BandContrastReader.o objs/Facets.o objs/FacetsWriter.o objs/FacetsReader.o objs/AFM.o objs/AFMDataReader.o objs/AFMDataWriter.o objs/Utility.o objs/Atom.o objs/SampleStretcher.o objs/Tilt.o objs/BitmapWriter.o objs/BandContrastAFMMApper.o objs/Amoeba.o

DEPENDENCIES = Vector3.o Matrix.o Quaternion.o BandContrast.o BandContrastWriter.o BandContrastReader.o Facets.o FacetsWriter.o FacetsReader.o AFM.o AFMDataReader.o AFMDataWriter.o Utility.o Atom.o SampleStretcher.o Tilt.o BitmapWriter.o BandContrastAFMMApper.o Amoeba.o

BandContrastSim: $(DEPENDENCIES)
	gcc -g -lm -o BandContrastSim $(OBJS) BandContrastSim.c
	ln -f BandContrastSim ~/prog/

Amoeba.o: BandContrastAFMMApper.o objs
	gcc -c -o objs/Amoeba.o Amoeba.c

SampleStretcher.o: Tilt.o objs
	gcc -c -o objs/SampleStretcher.o SampleStretcher.c

Tilt.o: AFM.o BandContrast.o objs
	gcc -c -o objs/Tilt.o Tilt.c

BandContrastReader.o: BandContrast.o objs
	gcc -c -o objs/BandContrastReader.o BandContrastReader.c

BandContrastWriter.o: BandContrast.o objs
	gcc -c -o objs/BandContrastWriter.o BandContrastWriter.c

BandContrast.o: Vector3.o Facets.o objs
	gcc -c -o objs/BandContrast.o BandContrast.c

BandContrastAFMMApper.o: BandContrast.o AFM.o objs
	gcc -c -o objs/BandContrastAFMMApper.o BandContrastAFMMApper.c
	
AFMDataWriter.o: AFM.o objs
	gcc -c -o objs/AFMDataWriter.o AFMDataWriter.c

AFMDataReader.o: AFM.o Utility.o objs
	gcc -c -o objs/AFMDataReader.o AFMDataReader.c

AFM.o: Utility.o objs
	gcc -c -o objs/AFM.o AFM.c

FacetsReader.o: Facets.o objs
	gcc -c -o objs/FacetsReader.o FacetsReader.c

FacetsWriter.o: Facets.o objs
	gcc -c -o objs/FacetsWriter.o FacetsWriter.c

Facets.o: AFM.o Matrix.o Vector3.o objs
	gcc -c -o objs/Facets.o Facets.c

Atom.o: Vector3.o Quaternion.o objs
	gcc -c -o objs/Atom.o Atom.c

Matrix.o: Utility.o objs
	gcc -c -o objs/Matrix.o Matrix.c

Utility.o: Vector3.o objs
	gcc -c -o objs/Utility.o Utility.c

Quaternion.o: Vector3.o objs
	gcc -c -o objs/Quaternion.o Quaternion.c

Vector3.o: objs
	gcc -c -o objs/Vector3.o Vector3.c

BitmapWriter.o: objs
	gcc -c -o objs/BitmapWriter.o BitmapWriter.c

objs:
	mkdir objs

zip: clean
	rm -fr ../BandContrastSim.zip
	zip -r ../BandContrastSim.zip .

clean: objs
	rm -r objs
	rm -f *.exe

test: Vector3.o
	gcc -g -lm -o test objs/Vector3.o test.c

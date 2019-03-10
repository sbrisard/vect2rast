all:
	gcc -c -std=c11 -Wall -Werror -fpic `pkg-config --cflags --libs glib-2.0` voxelize.c
	gcc -shared  voxelize.o `pkg-config --cflags --libs glib-2.0` -o libvoxelize.dll
	gcc -O2 test_voxelize.c -L. -lvoxelize `pkg-config --cflags --libs glib-2.0` -o test_voxelize.exe

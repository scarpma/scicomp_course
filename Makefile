## specifical rule hand-made
#trapezoidal_c.exe: trapezoidal.c
#	gcc trapezoidal.c -o trapezoidal_c.exe

# define the first target: all
# to build all, build all these targets
all: trapezoidal_c.exe trapezoidal_f.exe \
     trapezoidal_static_arr_c.exe trapezoidal_static_arr_f.exe \
     hello_c.exe

# generic rules for compiling .c and .F03 files
%_f.exe: %.F03
	gfortran $< -o $@

%_c.exe: %.c
	gcc -lm $< -o $@

clean:
	rm *.exe

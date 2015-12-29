CFLAGS  = -O3 -ffast-math   
CCFLAGS = -O3 -ffast-math -Wno-deprecated 
LIBS = -lm 

OBJ =  random.o pairlist.o mesofluid.o main.o
INC = *.h



#LDFLAGS = -s
C = gcc 
CC = g++

# Link rule

mydpd:	$(OBJ)
	$(CC) $(LDFLAGS) $(CCFLAGS) $(OBJ) $(LIBS) -o $@

clean:
	rm -f -r *.o mydpd 

# Compilation rules


.SUFFIXES: .o .cpp .C
.cpp.o:
	$(CC) $(CCFLAGS) -c $<
.C.o:
	$(CC) $(CCFLAGS) -c $<
.c.o:
	$(C) $(CFLAGS) -c $<


# Individual dependencies

$(OBJ):	$(INC)

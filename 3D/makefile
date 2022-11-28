XX=g++
SRCS=code/analysis.cpp\
     code/Node.cpp\
     code/Superball.cpp\
		 code/Vector.cpp\
		 code/Voxel.cpp\

OBJS=$(SRCS:.cpp=.o)

EXEC=analysis

start:$(OBJS)
	$(XX) -O3 -o $(EXEC) $(OBJS)

.cpp.o:
	$(XX) -O3 -o $@ -c $< 

clean:
	rm $(OBJS)
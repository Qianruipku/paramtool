CC=icc
BIN_PATH=/media/qianrui/work/LQR/code/bin
TEST=OFF
#-----------------------------------------------
OBJ=obj/main.o\
	obj/convert.o\
	obj/function.o\
	obj/mdpara.o\
	obj/elements.o\
	obj/conductivity.o\
	obj/chemical_potential.o\
	obj/tool.o

ifeq ($(TEST), ON)
    OPTION += -fsanitize=address -fno-omit-frame-pointer
    CC = g++
else
    OPTION += -D__COLOR -O3
endif

tool.exe:$(OBJ)
	$(CC) $(OPTION) -o tool.exe $(OBJ)
	@if [ $(TEST) == OFF ]; then cp tool.exe $(BIN_PATH); fi
	@if [ $(TEST) == ON ]; then cd test;./Autotest.sh ;cd ..; fi
obj/%.o:%.cpp
	$(CC) $(OPTION) -c $< -o $@
clean:
	@ rm -rf obj/*
	@ if [ ! -d tool.exe ]; then rm -rf tool.exe; fi 

CC=icc
BIN_PATH=/media/qianrui/work/LQR/code/bin
OBJ=obj/main.o\
	obj/convert.o\
	obj/function.o\
	obj/mdpara.o\
	obj/elements.o\
	obj/tool.o

tool.exe:$(OBJ)
	$(CC) -o tool.exe $(OBJ)
	mv tool.exe $(BIN_PATH)
obj/%.o:%.cpp
	$(CC) -c $< -o $@
clean:
	rm -rf obj/*

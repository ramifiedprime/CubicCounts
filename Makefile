#!/bin/bash
default:
	@echo "Compiling c_code/main.cn into main.exe..."
	@gcc -Wall c_code/main.c c_code/CCFCCFPRP.c -lm -o main.exe
	@echo "Complete."
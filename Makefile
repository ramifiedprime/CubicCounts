#!/bin/bash
default:
	@echo "Compiling c_code/main.cn into main.exe..."
	@gcc -Wall c_code/main.c c_code/CCFCCFPRP.c -lm -o CCFCCFPRP.exe
	@echo "Complete."

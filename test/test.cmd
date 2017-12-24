@echo off
	main.exe || exit
	stupid.exe <input.txt >output2.txt || exit
	check.exe || exit
	echo OK
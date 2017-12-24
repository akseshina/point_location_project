@echo off
	main.exe <test.in >test.out || exit
	stupid.exe <test.in >test2.out || exit
	check.exe || exit
	echo Test %%i - OK
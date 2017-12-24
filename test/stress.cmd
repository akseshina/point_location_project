@echo off

for /L %%i in (1,1,100000) do (
	gen.exe >input.txt || exit
	main.exe || exit
	stupid.exe <input.txt >output2.txt || exit
	check.exe || exit
	echo Test %%i - OK
)
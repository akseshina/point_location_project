@echo off

for /L %%i in (1,1,100000) do (
	gen.exe >test.in || exit
	echo test generated
	main.exe <test.in >test.out || exit
	echo main ended
	stupid.exe <test.in >test2.out || exit
	echo stupid ended
	check.exe || exit
	echo Test %%i - OK
)
all:
	@cmake -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES
clean:
	rm -rf bin build
full:
	@cmake -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -DENABLE_TESTS=ON
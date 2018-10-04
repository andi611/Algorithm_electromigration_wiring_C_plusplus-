# Algorithm_electromigration_wiring_C_plusplus
## The Wiring Topology for Electromigration Avoidance Problem: Optimal Wiring Topology Implementation in C++

## Description
- 使用之程式語言：**< C++ >**
- 使用之編譯器：**< g++ [gcc version 5.4.0 (GCC)] >**
- 各檔案說明：
	- main.cpp : C++ source code 主程式
	- pa2  : main.cpp compile 後的 Unix executable binary code執行檔
	- readme.md : 本說明文字檔
	- .in/.out : 測資input及output檔

- 編譯方式說明：        	
    * 主程式：main.cpp
		- 請在主程式的目錄下，鍵入以下指令，即可完成編譯：
		- g++ -O3 main.cpp -std=c++11
		- 在主程式的目錄下會產生一個名為 a.out 的執行檔
		- 如果要命名編譯後的執行檔，請先鍵入以下指令：
		- g++ -O3 main.cpp -std=c++11 -o [my_name]
		- 如：g++ -O3 main.cpp -std=c++11 -o pa2

- 執行、使用方式說明：
   	* 主程式：
    	- 編譯完成後，在檔案目錄下會產生一個 a.out 的執行檔
   	* 執行檔的命令格式為：
	   	- ./a.out <input file name> <output file name>
	   	- 例如：./a.out test*.in test*.out
	   	- 若是命名過後的執行檔(如pa2)，則執行檔的命令格式為：
		- ./my_name <input file name> <output file name>
	 	- 例如：./pa2 test*.in test*.out
     
- 執行結果說明（說明執行結果的觀看方法，及解釋各項數據等）：
   	* 主程式：
		- 主程式執行時會在 console 輸出讀到的檔案資訊，以及目前執行的步驟。
		- 例如：
		- 執行 ./pa2 INP1.txt 1.out時，console會輸出以下內容：
		
		- Number of sources + sinks: 7
		- Number of source: 3
		- Number of sink: 4
		- Complete building flow network with total flow: 19
		- Pushing initial flow with normal Method.
		- Building residual network, starting negative cycle removal...
		- Iteration: 1
		- Final Area: 142
		- Result successfully saved to: 1.out
